/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2022 Mathieu Fortin for Rouge-Epicea
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3 of the License, or (at your option) any later version.
 *
 * This library is distributed with the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * Please see the license at http://www.gnu.org/copyleft/lesser.html.
 */
package repicea.stats.model.glm.measerr;

import java.security.InvalidParameterException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import repicea.math.Matrix;
import repicea.stats.Distribution;
import repicea.stats.StatisticalUtility;
import repicea.stats.data.DataSet;
import repicea.stats.data.GenericStatisticalDataStructure;
import repicea.stats.data.StatisticalDataException;
import repicea.stats.estimates.Estimate;
import repicea.stats.estimators.AbstractEstimator.EstimatorCompatibleModel;
import repicea.stats.estimators.Estimator;
import repicea.stats.model.AbstractStatisticalModel;
import repicea.stats.model.CompositeLogLikelihoodWithExplanatoryVariables;
import repicea.stats.model.IndividualLogLikelihood;
import repicea.stats.model.PredictableModel;
import repicea.stats.model.glm.Family.GLMDistribution;
import repicea.stats.model.glm.GeneralizedLinearModel;
import repicea.stats.model.glm.LinkFunction.Type;

/**
 * This class implements the SIMEX method.
 * @author Mathieu Fortin - August 2022
 */
public class SIMEXModel extends AbstractStatisticalModel implements EstimatorCompatibleModel, PredictableModel {

	/**
	 * A fake GLM just to produce the predictions relying on the SIMEX parameter estimates.
	 * @author Mathieu Fortin - September 2022
	 */
	final class PredGLM extends GeneralizedLinearModel {

		/**
		 * A fake estimator for the generalized linear model that is used to 
		 * generate the predictions.
		 * @author Mathieu Fortin - September 2022
		 */
		final class FakeEstimator implements Estimator {

			Estimate<?> parameterEstimates;
			
			@Override
			public boolean doEstimation() throws EstimatorException {return true;}

			@Override
			public boolean isConvergenceAchieved() {return true;}

			@Override
			public Estimate<?> getParameterEstimates() {return parameterEstimates;}

			@Override
			public DataSet getParameterEstimatesReport() {return null;}

			@Override
			public DataSet getConvergenceStatusReport() {return null;}
			
		}
		
		protected PredGLM(GeneralizedLinearModel glm) {
			super(glm);
			this.setEstimator(new FakeEstimator());
		}

		@Override
		public FakeEstimator getEstimator() {
			return (FakeEstimator) super.getEstimator();
		}

	}
	
	
	final class InternalGLM extends GeneralizedLinearModel implements Cloneable {
		
		private final class SIMEXDataStructure extends GenericStatisticalDataStructure {

			Matrix additionalMeasErr;
			double sqrtFactor;
			
			private SIMEXDataStructure(DataSet dataSet) {
				super(dataSet);
				additionalMeasErr = new Matrix(dataSet.getNumberOfObservations(), 1);
				sqrtFactor = 0d;
			}
			
			@Override
			protected Matrix getVectorOfThisField(String fName) {
				Matrix originalValue = super.getVectorOfThisField(fName);
				return varWithMeasErr.equals(fName) ?  // we only add the random deviate if this is the field with measurement error
						originalValue.add(additionalMeasErr.scalarMultiply(sqrtFactor)) :
							originalValue;
			}
		}

		/**
		 * Main constructor.
		 * @param glm
		 */
		private InternalGLM(GeneralizedLinearModel glm) {
			super(glm);
		}

		/**
		 * Constructor for clones.
		 * @param dataSet
		 * @param linkFunctionType
		 * @param modelDefinition
		 * @param startingParms
		 */
		private InternalGLM(DataSet dataSet, GLMDistribution d, Type linkFunctionType, String modelDefinition, Matrix startingParms) {
			super(dataSet, d, linkFunctionType, modelDefinition);
		}

		@Override
		protected SIMEXDataStructure createDataStructure(DataSet dataSet, Object addParm) {
			return new SIMEXDataStructure(dataSet);
		}

		@Override
		protected SIMEXDataStructure getDataStructure() {
			return (SIMEXDataStructure) super.getDataStructure();
		}

		@Override
		public InternalGLM clone() {
			InternalGLM clone = new InternalGLM(this.getDataStructure().getDataSet(), 
					this.getDistribution(),
					this.getLinkFunctionType(), 
					this.getModelDefinition(), null);
			clone.getCompleteLogLikelihood().variance = this.getCompleteLogLikelihood().variance.getDeepClone();
			return clone;
		}
		
		@Override
		protected BootstrapCompositeLogLikelihoodWithExplanatoryVariables createCompleteLLK(Object addParm) {
			return new BootstrapCompositeLogLikelihoodWithExplanatoryVariables(individualLLK, y, this);
		}

		@Override
		public BootstrapCompositeLogLikelihoodWithExplanatoryVariables getCompleteLogLikelihood() {
			return (BootstrapCompositeLogLikelihoodWithExplanatoryVariables) super.getCompleteLogLikelihood();
		}
	}
	
	
	@SuppressWarnings("serial")
	static final class BootstrapCompositeLogLikelihoodWithExplanatoryVariables extends CompositeLogLikelihoodWithExplanatoryVariables {

		final InternalGLM caller;
		Matrix variance;
		Matrix std;
//		Matrix xTmp;
//		int indexVarWithMeasErr;
		
		private BootstrapCompositeLogLikelihoodWithExplanatoryVariables(IndividualLogLikelihood innerLogLikelihoodFunction, Matrix yValues, InternalGLM caller) {
			super(innerLogLikelihoodFunction, null, yValues);
			this.caller = caller;
		}
		
		void generateMeasurementError(double factor) {
			this.caller.getDataStructure().sqrtFactor = Math.sqrt(factor);
			if (std == null) {
				std = variance.elementWisePower(0.5);
			}
			if (OverrideVarianceForTest) {
				this.caller.getDataStructure().additionalMeasErr = StatisticalUtility.drawRandomVector(std.m_iRows, Distribution.Type.GAUSSIAN).scalarMultiply(1.59);
			} else {
				this.caller.getDataStructure().additionalMeasErr = std.elementWiseMultiply(StatisticalUtility.drawRandomVector(std.m_iRows, Distribution.Type.GAUSSIAN));
			}

			xValues = this.caller.getDataStructure().constructMatrixX();
//			xTmp = this.xValues.getDeepClone();
//			Matrix additionalMeasErr;
//			for (int i = 0; i < xTmp.m_iRows; i++) {
//				double currentValue = xTmp.getValueAt(i, indexVarWithMeasErr);
//				xTmp.setValueAt(i, indexVarWithMeasErr, currentValue + additionalMeasErr.getValueAt(i, 0) * sqrtFactor);
//			}
		}
		
//		protected void setValuesInLikelihoodFunction(int index) {
//			super.setValuesInLikelihoodFunction(index);
//			getOriginalFunction().setVariables(xTmp.getSubMatrix(index, index, 0, xTmp.m_iCols - 1));
//		}
 
		
	}
	
	static boolean OverrideVarianceForTest = false;
	
	protected final InternalGLM originalGLM;
	protected final PredGLM predGLM;
	protected final String varWithMeasErr;
	int nbBootstrapRealizations = 100;
	int nbThreads = 2;
	double[] factors = new double[] {0, .2, .4, .6, .8, 1, 1.2, 1.4, 1.6, 1.8, 2.0};
	
	
	
	/**
	 * Constructor.
	 * @param glm a GeneralizedLinearModel instance
	 * @param varWithMeasErr the field name of the variable with measurement error
	 * @param varianceField the field name of the variance of the measurement error
	 */
	public SIMEXModel(GeneralizedLinearModel glm, String varWithMeasErr, String varianceField) {
		this.varWithMeasErr = varWithMeasErr;
		this.originalGLM = new InternalGLM(glm);
		this.predGLM = new PredGLM(glm);
		try {
			setModelDefinition(glm.getModelDefinition());
		} catch (StatisticalDataException e) {
			e.printStackTrace();
		}
		int indexVarWithMeasErr = originalGLM.getDataStructure().indexOfThisEffect(varWithMeasErr);
		if (indexVarWithMeasErr == -1)
			throw new InvalidParameterException("The effect " + varWithMeasErr + " is not part of the model definition!");
		DataSet ds = originalGLM.getDataStructure().getDataSet();
		int indexVarianceField = ds.getIndexOfThisField(varianceField);
		if (indexVarianceField == -1)
			throw new InvalidParameterException("The field " + varianceField + " is not found in the original dataset!");
		Matrix varVector = new Matrix(ds.getNumberOfObservations(), 1);
		for (int i = 0; i < ds.getNumberOfObservations(); i++) {
			varVector.setValueAt(i, 0, ((Number) ds.getObservations().get(i).getValueAt(indexVarianceField)).doubleValue());
		}
		this.originalGLM.getCompleteLogLikelihood().variance = varVector;
	}
	
	/**
	 * Set the number of bootstrap realizations. <br>
	 * <br>
	 * By default, the number of realizations is set to 100.
	 * @param nbReal an integer between 1 and 1000.
	 */
	public void setNumberOfBootstrapRealizations(int nbReal) {
		if (nbReal < 1 || nbReal > 1000) {
			throw new InvalidParameterException("The number of bootstrap realizations must range between 1 and 1000!");
		}
		this.nbBootstrapRealizations = nbReal;
	}

	/**
	 * Provide the number of bootstrap realizations to be used in the simulation step.
	 * @return an integer
	 */
	public int getNumberOfBootstrapRealizations() {
		return this.nbBootstrapRealizations;
	}
	
	/**
	 * Set the number of threads. <br>
	 * <br>
	 * By default, the number of threads is set to 2.
	 * @param nbThreads an integer between 1 and 10.
	 */
	public void setNbThreads(int nbThreads) {
		if (nbThreads < 1 || nbThreads > 10) {
			throw new InvalidParameterException("The number of threads must range between 1 and 10!");
		}
		this.nbThreads = nbThreads;
	}

	/**
	 * Return the number of threads to be used in the simulation.
	 * @return an integer
	 */
	public int getNbThreads() {return nbThreads;}
	
	/**
	 * Set the variance inflation factors for the simulation. <br>
	 * <br>
	 * By default the factors are set to 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0 .
	 * 
	 * @param factors a series of doubles ranging from 0 to 2.
	 */
	public void setFactors(double...factors) {
		if (factors == null || factors.length == 0) {
			throw new InvalidParameterException("The factors argument must be non null!");
		}
		for (double d : factors) {
			if (d < 0 || d > 2) {
				throw new InvalidParameterException("The factors must range between 0 and 2!");
			}
		}
		this.factors = factors;
	}
	
	/**
	 * Return a copy of the current factors.<br>
	 * <br>
	 * Changing the resulting array has no impact since it is a copy of the original. To change the
	 * factors, use the setFactors() method instead.
	 * 
	 * @return an array of doubles
	 */
	public double[] getFactors() {
		return Arrays.copyOf(factors, factors.length);
	}
	
	@Override
	public List<String> getOtherParameterNames() {
		List<String> names = new ArrayList<String>();
		names.addAll(this.originalGLM.getOtherParameterNames());
		return names;
	}
	
	@Override
	public void setEstimator(Estimator e) {
		if (e instanceof SIMEXEstimator) {
			this.estimator = e;
		} else {
			throw new UnsupportedOperationException("The SIMEX model does not allow for other estimators than its default estimator!");
		} 
	}
	
	@Override
	public String toString() {return "SIMEX model";}
	

	@Override
	protected Estimator instantiateDefaultEstimator() {return new SIMEXEstimator(this);}

	@Override
	public boolean isInterceptModel() {return originalGLM.isInterceptModel();}

	@Override
	public List<String> getEffectList() {return originalGLM.getEffectList();}

	@Override
	public int getNumberOfObservations() {return originalGLM.getNumberOfObservations();}

	/**
	 * Provide a dataset with the observed parameter estimates for each level of inflated variance. <br>
	 * <br>
	 * This dataset is only available if the estimator has converged. If not, the method returns null.
	 * @return a DataSet instance (or null if the estimator has not converged)
	 */
	public DataSet getObservedParameterEstimates() {
		return getEstimator().isConvergenceAchieved() ?
			((SIMEXEstimator) getEstimator()).parmsObsDataSet :
				null;
	}

	/**
	 * Provide a dataset with the predicted parameter estimates for each level of inflated variance. <br>
	 * <br>
	 * The predicted parameter estimates result from the regression of the observed parameter estimates
	 * against the level of inflated variance. The predictions are generated for the range -1 to the maximum
	 * of the factors, by 0.01. <br>
	 * <br>
	 * This dataset is only available if the estimator has converged. If not, the method returns null.
	 * @return a DataSet instance (or null if the estimator has not converged)
	 */
	public DataSet getPredictedParameterEstimates() {
		return getEstimator().isConvergenceAchieved() ?
			((SIMEXEstimator) getEstimator()).parmsPredDataSet :
				null;
	}

	@Override
	public Matrix getPredicted() {return predGLM.getPredicted();}

	@Override
	public Matrix getResiduals() {return predGLM.getResiduals();}

}
