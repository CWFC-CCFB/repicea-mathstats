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
import java.util.List;

import repicea.math.AbstractMathematicalFunction;
import repicea.math.AbstractMathematicalFunctionWrapper;
import repicea.math.InternalMathematicalFunctionWrapper;
import repicea.math.MathematicalFunction;
import repicea.math.Matrix;
import repicea.math.ProductFunctionWrapper;
import repicea.math.SymmetricMatrix;
import repicea.stats.data.DataSet;
import repicea.stats.data.GenericStatisticalDataStructure;
import repicea.stats.data.StatisticalDataStructure;
import repicea.stats.distributions.utility.GaussianUtility;
import repicea.stats.functions.EmpiricalDistributionProbabilityDensityFunction;
import repicea.stats.integral.TrapezoidalRule;
import repicea.stats.model.CompositeLogLikelihoodWithExplanatoryVariables;
import repicea.stats.model.IndividualLogLikelihood;
import repicea.stats.model.Likelihood;
import repicea.stats.model.WrappedIndividualLogLikelihood;
import repicea.stats.model.glm.BernoulliIndividualLikelihood;
import repicea.stats.model.glm.LinkFunction;
import repicea.stats.model.glm.LinkFunction.Type;

public class GLMNormalClassicalMeasErrorDefinition extends AbstractGLMMeasErrorDefinition {

	private static class GLMWithNormalClassicalMeasErrorDataStructure extends GenericStatisticalDataStructure {

		private final GLMNormalClassicalMeasErrorDefinition measErr;
		
		private GLMWithNormalClassicalMeasErrorDataStructure(DataSet dataSet, GLMNormalClassicalMeasErrorDefinition measErr) {
			super(dataSet);
			this.measErr = measErr;
		}
		
		void initializeMeasurementErrorDefinition() {
			measErr.wVector = new ArrayList<Double>(); 
			// populating the wVector
			int indexField = dataSet.getIndexOfThisField(measErr.effectWithMeasError);
			for (int i = 0; i < dataSet.getNumberOfObservations(); i++) {
				measErr.wVector.add((Double) dataSet.getObservations().get(i).getValueAt(indexField));
			}
			measErr.varianceVector = new ArrayList<Double>();
			// populating variance
			indexField = dataSet.getIndexOfThisField(measErr.varianceField);
			for (int i = 0; i < dataSet.getNumberOfObservations(); i++) {
				measErr.varianceVector.add((Double) dataSet.getObservations().get(i).getValueAt(indexField));
			}
		}

	}
	
	@SuppressWarnings("serial")
	static class AlteredGaussianFunction extends AbstractMathematicalFunction {

		static final int SIGMA2_INDEX = 1;
		static final int X_INDEX = 0;

		double wValue;
		
		AlteredGaussianFunction() {
			setVariableValue(X_INDEX, 0d);
			setVariableValue(SIGMA2_INDEX, 1d);
		}

		@Override
		public void setParameterValue(int index, double value) {
			throw new UnsupportedOperationException("The AlteredGaussianFunction only has variables!");
		}

		@Override
		public void setVariableValue(int index, double value) {
			if (index > 1) {
				throw new InvalidParameterException("The AlteredGaussianFunction only has two variable (namely the observation x and the variance)!");
			} else {
				if (index == SIGMA2_INDEX && value <= 0) {
					throw new InvalidParameterException("The sigma2 variable must be strictly positive!");
				}
				super.setVariableValue(index, value);
			}
		}
		
		@Override
		public Double getValue() {
			double mu = getVariableValue(X_INDEX);
			double sigma2 = getVariableValue(SIGMA2_INDEX);
			return GaussianUtility.getProbabilityDensity(wValue, mu, sigma2);
		}

		@Override
		public Matrix getGradient() {
			return null;
		}

		@Override
		public SymmetricMatrix getHessian() {
			return null;
		}

	}

	
//	@SuppressWarnings("serial")
//	static class InternalProductFunctionWrapper extends ProductFunctionWrapper implements GaussHermiteQuadratureCompatibleFunction<Double> {
//
//		private final AlteredGaussianFunction agf;
//		
//		InternalProductFunctionWrapper(AlteredGaussianFunction agf, InternalMathematicalFunctionWrapper... wrappedOriginalFunctions) {
//			super(wrappedOriginalFunctions);
//			this.agf = agf;
//		}
//		
//		@Override
//		public double convertFromGaussToOriginal(double x, double mu, int covarianceIndexI, int covarianceIndexJ) {
//			return agf.wValue - Math.sqrt(2 * agf.getParameterValue(AlteredGaussianFunction.SIGMA2_INDEX)) * x;
//		}
//
//		@Override
//		public double getIntegralAdjustment(int dimensions) {
//			return - Math.sqrt(2 * agf.getParameterValue(AlteredGaussianFunction.SIGMA2_INDEX));
//		}
//	}
	
	@SuppressWarnings("serial")
	static class LikelihoodGLMWithNormalClassicalMeasErr extends AbstractMathematicalFunctionWrapper implements Likelihood {

		private final GLMNormalClassicalMeasErrorDefinition measErr;
		private final TrapezoidalRule integrator;
		private final GradientHessianProvider gradientProvider;
		private final GradientHessianProvider hessianProvider;


		public LikelihoodGLMWithNormalClassicalMeasErr(ProductFunctionWrapper pfw, GLMNormalClassicalMeasErrorDefinition measErr) {
			super(pfw);
			this.measErr = measErr;
			integrator = new TrapezoidalRule(measErr.resolution);
			gradientProvider = new GradientHessianProvider(pfw , true) ;
			hessianProvider = new GradientHessianProvider(pfw, false);
		}

		@Override
		public void setYVector(Matrix yVector) {
			measErr.originalModelLikelihood.setYVector(yVector);
		}

		@Override
		public Matrix getYVector() {return null;}

		@Override
		public Matrix getPredictionVector() {return null;}

		@Override
		public ProductFunctionWrapper getOriginalFunction() {
			return (ProductFunctionWrapper) super.getOriginalFunction();
		}
		
		@Override
		public Double getValue() {
			double integralApprox = integrator.getIntegralApproximation(getOriginalFunction(), measErr.indexEffectWithMeasError, false);
			return integralApprox;
		}

		@Override
		public Matrix getGradient() {
			return integrator.getIntegralApproximationForMatrixFunction(gradientProvider, measErr.indexEffectWithMeasError, false);
		}

		@Override
		public SymmetricMatrix getHessian() {
			SymmetricMatrix hessian = SymmetricMatrix.convertToSymmetricIfPossible(integrator.getIntegralApproximationForMatrixFunction(hessianProvider, measErr.indexEffectWithMeasError, false));
			return hessian;
		}

	}
	

	
	@SuppressWarnings("serial")
	static class GLMWithNormalClassicalMeasErrorCompositeLogLikelihood extends CompositeLogLikelihoodWithExplanatoryVariables {
		
		private final GLMNormalClassicalMeasErrorDefinition measErr;
		
		public GLMWithNormalClassicalMeasErrorCompositeLogLikelihood(IndividualLogLikelihood innerLogLikelihoodFunction, 
				Matrix xValues, 
				Matrix yValues, 
				GLMNormalClassicalMeasErrorDefinition measErr) {
			super(innerLogLikelihoodFunction, xValues, yValues);
			this.measErr = measErr;
		}
		
		@Override
		protected void setValuesInLikelihoodFunction(int index) {
			measErr.setValuesForObservation(index);
			super.setValuesInLikelihoodFunction(index);
		}
	}

	
	private final double resolution;
	private List<Double> wVector;
	private List<Double> varianceVector;
	
	private BernoulliIndividualLikelihood originalModelLikelihood;
	private AlteredGaussianFunction errorModelLikelihood;
	private MathematicalFunction xDistribution;
	private final String varianceField;
	private LikelihoodGLMWithNormalClassicalMeasErr lk;
	private double maxX;
	
	public GLMNormalClassicalMeasErrorDefinition(String effectWithMeasError, String varianceField, double resolution) {
		super(MeasurementErrorModel.Classical, effectWithMeasError);
		if (resolution <= 0) {
			throw new InvalidParameterException("The resolution argument must be strictly positive!");
		} else {
			this.resolution = resolution;
		}
		this.varianceField = varianceField;
	}

	private void setValuesForObservation(int index) {
		errorModelLikelihood.wValue = wVector.get(index);
		double sigma2 = varianceVector.get(index);
		errorModelLikelihood.setVariableValue(AlteredGaussianFunction.SIGMA2_INDEX, sigma2);
		double p = 1E-16;
		double delta = Math.sqrt(-2 * sigma2 * Math.log(Math.sqrt(2 * Math.PI * sigma2) * p));
		double lowerBound = errorModelLikelihood.wValue - delta;
		double upperBound = errorModelLikelihood.wValue + delta;
		if (xDistribution instanceof EmpiricalDistributionProbabilityDensityFunction) {
			double officialLowerBound = ((EmpiricalDistributionProbabilityDensityFunction) xDistribution).getLowerBoundForX();
			if (lowerBound <= officialLowerBound) 
				lowerBound = officialLowerBound;
			double officialUpperBound = ((EmpiricalDistributionProbabilityDensityFunction) xDistribution).getUpperBoundForX();
			if (upperBound > officialUpperBound) 
				upperBound = officialUpperBound;
		} else {
			if (lowerBound <= 0) 
				lowerBound = 0.1;
			if (upperBound > maxX) 
				upperBound = maxX;
		}
		lk.integrator.setLowerBound(lowerBound);
		lk.integrator.setUpperBound(upperBound);
	}

	/**
	 * Approximate the first two central moments of the distribution of X from the sample of w. <br>
	 * <br>
	 * These values are used as starting values.
	 * 
	 * @return an 2-slot array (the first is the mean, the second is the variance)
	 */
	private double[] getLogScaleMuAndVarianceAndMaximumX() {
		double[] muAndSigma2 = new double[3];
		double maxX = Double.NEGATIVE_INFINITY;
		double sumW = 0d;
		for (Double w : wVector) {
			sumW += Math.log(w);
			maxX = w > maxX ? w : maxX;
		}
		muAndSigma2[0] = sumW / wVector.size();
		muAndSigma2[2] = maxX;
		double sse = 0d;
		for (Double w : wVector) {
			sse += (Math.log(w) - muAndSigma2[0]) * (Math.log(w) - muAndSigma2[0]);
		}
		muAndSigma2[1] = sse / (wVector.size() - 1);
		return muAndSigma2;
	}
	
	@Override
	public void validate(GLMWithMeasurementError glm) {
		super.validate(glm);
		((GLMWithNormalClassicalMeasErrorDataStructure) glm.getDataStructure()).initializeMeasurementErrorDefinition();
	}

	@Override
	public CompositeLogLikelihoodWithExplanatoryVariables createCompositeLikelihoodFromModel(GLMWithMeasurementError glm) {
		return new GLMWithNormalClassicalMeasErrorCompositeLogLikelihood(glm.getIndividualLogLikelihood(), 
				glm.getMatrixX(), 
				glm.getVectorY(), 
				this);
	}

	@Override
	public StatisticalDataStructure createDataStructureFromDataSet(DataSet dataSet) {
		return new GLMWithNormalClassicalMeasErrorDataStructure(dataSet, this);
	}

	@Override
	public IndividualLogLikelihood createIndividualLogLikelihoodFromModel(GLMWithMeasurementError glm) {
		originalModelLikelihood = new BernoulliIndividualLikelihood(glm.getLinkFunction());
		InternalMathematicalFunctionWrapper wrapper1 = new InternalMathematicalFunctionWrapper(originalModelLikelihood, 
				InternalMathematicalFunctionWrapper.produceListFromTo(0, originalModelLikelihood.getNumberOfParameters() - 1),
				InternalMathematicalFunctionWrapper.produceListFromTo(0, originalModelLikelihood.getNumberOfVariables() - 1)); 
		errorModelLikelihood = new AlteredGaussianFunction();
		List<Integer> newVariableIndices = new ArrayList<Integer>(); 
		newVariableIndices.add(indexEffectWithMeasError);
		newVariableIndices.add(originalModelLikelihood.getNumberOfVariables());
		InternalMathematicalFunctionWrapper wrapper2 = new InternalMathematicalFunctionWrapper(errorModelLikelihood, 
				new ArrayList<Integer>(),
				newVariableIndices);
//		double[] approxMuAndSigma2 = getLogScaleMuAndVarianceAndMaximumX();
//		maxX = approxMuAndSigma2[2];
//		xDistribution = new LogGaussianFunction(approxMuAndSigma2[0], approxMuAndSigma2[1]);
//		InternalMathematicalFunctionWrapper wrapper3 = new InternalMathematicalFunctionWrapper(xDistribution, 
//				InternalMathematicalFunctionWrapper.produceListFromTo(originalModelLikelihood.getNumberOfParameters() + errorModelLikelihood.getNumberOfParameters(), 
//						originalModelLikelihood.getNumberOfParameters()  + errorModelLikelihood.getNumberOfParameters() + xDistribution.getNumberOfParameters() - 1),
//				InternalMathematicalFunctionWrapper.produceListFromTo(indexEffectWithMeasError, indexEffectWithMeasError)); 
		xDistribution = new EmpiricalDistributionProbabilityDensityFunction(wVector, null);	// null: no weighting here
		InternalMathematicalFunctionWrapper wrapper3 = new InternalMathematicalFunctionWrapper(xDistribution, 
				new ArrayList<Integer>(),
				InternalMathematicalFunctionWrapper.produceListFromTo(indexEffectWithMeasError, indexEffectWithMeasError)); 
		ProductFunctionWrapper pfw = new ProductFunctionWrapper(wrapper1, wrapper2, wrapper3);
		lk = new LikelihoodGLMWithNormalClassicalMeasErr(pfw, this);
		return new WrappedIndividualLogLikelihood(lk);
	}

	@Override
	public LinkFunction createLinkFunction(Type linkFunctionType, GLMWithMeasurementError glm) {
		LinkFunction lf = new LinkFunction(linkFunctionType);
		int dimension = glm.getMatrixX().m_iCols;
		lf.setParameters(new Matrix(dimension, 1));
		lf.setVariables(new Matrix(1, dimension));
		return lf;
	}

	@Override
	public List<String> getAdditionalEffects() {
		List<String> additionalEffects = new ArrayList<String>();
		additionalEffects.add("Variance error model");
		additionalEffects.add("Mean distribution of X");
		additionalEffects.add("Variance distribution of X");
		return additionalEffects;
	}
	

}
