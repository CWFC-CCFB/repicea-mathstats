/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2012 Mathieu Fortin for Rouge-Epicea
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
package repicea.stats.model.glm;


import java.security.InvalidParameterException;
import java.util.ArrayList;
import java.util.List;

import repicea.math.MathematicalFunction;
import repicea.math.Matrix;
import repicea.stats.data.DataSet;
import repicea.stats.data.GenericStatisticalDataStructure;
import repicea.stats.data.StatisticalDataException;
import repicea.stats.data.StatisticalDataStructure;
import repicea.stats.estimators.Estimator;
import repicea.stats.estimators.MaximumLikelihoodEstimator;
import repicea.stats.estimators.MaximumLikelihoodEstimator.MaximumLikelihoodCompatibleModel;
import repicea.stats.model.AbstractStatisticalModel;
import repicea.stats.model.CompositeLogLikelihoodWithExplanatoryVariables;
import repicea.stats.model.IndividualLikelihood;
import repicea.stats.model.IndividualLogLikelihood;
import repicea.stats.model.PredictableModel;
import repicea.stats.model.WrappedIndividualLogLikelihood;
import repicea.stats.model.glm.Family.GLMDistribution;
import repicea.stats.model.glm.LinkFunction.Type;

/**
 * This class implements generalized linear models. 
 * @author Mathieu Fortin - August 2011
 */
@SuppressWarnings("serial")
public class GeneralizedLinearModel extends AbstractStatisticalModel implements MaximumLikelihoodCompatibleModel, PredictableModel {

	static abstract class GLMIndividualLikelihood extends IndividualLikelihood {

		protected final List<Integer> additionalParameterIndices;

		protected GLMIndividualLikelihood(MathematicalFunction originalFunction) {
			super(originalFunction);
			additionalParameterIndices = new ArrayList<Integer>();
		}
		
		/**
		 * Record an index as one of an additional parameter. <br>
		 * <br>
		 * For instance, this method can be used to record the index of the dispersion 
		 * parameter in a negative binomial regression. 
		 * 
		 * @param index
		 */
		void recordAdditionalParameterIndex(int index) {
			additionalParameterIndices.add(index);
		}

	}
	
	private final StatisticalDataStructure dataStruct;
	private final CompositeLogLikelihoodWithExplanatoryVariables completeLLK;
	protected final IndividualLogLikelihood individualLLK;
//	protected final LinkFunction lf;
	protected final Family family;
	protected Matrix matrixX;		// reference
	protected Matrix y;				// reference
	
	private double convergenceCriterion = 1E-8;

	

	/**
	 * General constructor
	 * @param dataSet the fitting data
	 * @param d a Distribution enum that defines the distribution of the response variable
	 * @param linkFunctionType the type of link function (Logit, CLogLog, ...)
	 * @param modelDefinition a String that defines the dependent variable and the effects of the model
	 * @param startingBeta the starting values of the parameters
	 */
	public GeneralizedLinearModel(DataSet dataSet, GLMDistribution d, Type linkFunctionType, String modelDefinition, Matrix startingBeta) {
		this(dataSet, d, linkFunctionType, modelDefinition, null, startingBeta, null);
	}
	
	/**
	 * Generic private constructor.
	 * @param dataSet a DataSet instance
	 * @param d a GLMDistribution enum that defines the distribution of the response variable
	 * @param linkFunctionType a Type enum that defines the type of linkFunction
	 * @param modelDefinition a String that defines the model 
	 * @param llk an IndividualLogLikelihood instance 
	 * @param startingBeta a Matrix of starting parameters for the fixed effects
	 * @param additionalParm an object that contains additional parameters (can be null)
	 */
	protected GeneralizedLinearModel(DataSet dataSet, GLMDistribution d, Type linkFunctionType, String modelDefinition, IndividualLogLikelihood llk, Matrix startingBeta, Object additionalParm) {
		super();
		dataStruct = createDataStructure(dataSet, additionalParm);

		// then define the model effects and retrieve matrix X and vector y
		try {
			setModelDefinition(modelDefinition, additionalParm);
		} catch (StatisticalDataException e) {
			System.out.println("Unable to define this model : " + modelDefinition);
			e.printStackTrace();
		}
		family = createFamily(d, linkFunctionType, additionalParm);
		if (llk != null) {
			individualLLK = llk;
		} else {
			individualLLK = createIndividualLLK(additionalParm);
		}
		completeLLK = createCompleteLLK(additionalParm);
		Matrix startingParms = startingBeta == null ?
				d.getStartingParms(matrixX.m_iCols) :
					d.getStartingParms(startingBeta);
		for (int i = matrixX.m_iCols; i < startingParms.m_iRows; i++) {
			GLMIndividualLikelihood glmLk = ((GLMIndividualLikelihood) ((WrappedIndividualLogLikelihood) individualLLK).getOriginalFunction());
			glmLk.recordAdditionalParameterIndex(i);
		}
		setParameters(startingParms);
	}
	
	
	/**
	 * Constructor using a vector of 0s as starting values for the parameters
	 * @param dataSet the fitting data
	 * @param d a Distribution enum that defines the distribution of the response variable
	 * @param linkFunctionType the type of link function (Logit, CLogLog, ...)
	 * @param modelDefinition a String that defines the dependent variable and the effects of the model
	 */
	public GeneralizedLinearModel(DataSet dataSet, GLMDistribution d, Type linkFunctionType, String modelDefinition) {
		this(dataSet, d, linkFunctionType, modelDefinition, null);
	}

	/**
	 * Constructor for derived class.
	 * @param glm a GeneralizedLinearModel instance
	 */
	protected GeneralizedLinearModel(GeneralizedLinearModel glm) {
		this(glm.getDataStructure().getDataSet(), 
				glm.family.dist, 
				glm.family.lf.getType(), 
				glm.getModelDefinition(), 
				glm.individualLLK, 
				null, 
				null);
	}


	protected StatisticalDataStructure createDataStructure(DataSet dataSet, Object addParm) {
		return new GenericStatisticalDataStructure(dataSet);
	}

	@Override
	public List<String> getOtherParameterNames() {
		List<String> names = new ArrayList<String>();
		names.addAll(family.dist.additionalParmNames);
		return names;
	}

	/**
	 * Provide the data structure underlying this model instance.
	 * @return a StatisticalDataStructure instance
	 */
	protected StatisticalDataStructure getDataStructure() {
		return (GenericStatisticalDataStructure) dataStruct;
	}
	
	protected CompositeLogLikelihoodWithExplanatoryVariables createCompleteLLK(Object addParm) {
		return new CompositeLogLikelihoodWithExplanatoryVariables(individualLLK, matrixX, y);
	}

	protected Family createFamily(GLMDistribution d, Type linkFunctionType, Object addParm) {
		return Family.createFamily(d, linkFunctionType, null);
//		return new LinkFunction(linkFunctionType);
	}

	protected IndividualLogLikelihood createIndividualLLK(Object addParm) {
		GLMIndividualLikelihood indLk;
		switch(family.dist) {
		case Bernoulli:
			indLk = new BernoulliIndividualLikelihood(family.lf);
			break;
		case NegativeBinomial:
			indLk = new NegativeBinomialIndividualLikelihood(family.lf);
			break;
		default:
			throw new InvalidParameterException("The distribution " + family.dist.name() + " is not supported yet!"); 
		}
		return new WrappedIndividualLogLikelihood(indLk);
	}

	/**
	 * Provide the distribution of the response variable.
	 * @return a Distribution enum
	 */
	public GLMDistribution getDistribution() {return family.dist;}
	
	/**
	 * This method returns the type of the link function.
	 * @return a LinkFunction.Type enum variable
	 */
	public LinkFunction.Type getLinkFunctionType() {
		return getLinkFunction().getType();
	}
	
	protected LinkFunction getLinkFunction() {return family.lf;}

	@Override
	protected void setModelDefinition(String modelDefinition, Object additionalParm) throws StatisticalDataException {
		super.setModelDefinition(modelDefinition, additionalParm);
		getDataStructure().setModelDefinition(modelDefinition);
		matrixX = getDataStructure().constructMatrixX();
		y = getDataStructure().constructVectorY();
	}
	

	@Override
	protected Estimator instantiateDefaultEstimator() {
		return new MaximumLikelihoodEstimator(this);
	}

	@Override
	public String toString() {
		return "Generalized linear model";
	}

	@Override
	public int getNumberOfObservations() {return getDataStructure().getNumberOfObservations();}
	
	public Matrix getLinearPredictions() {
		if (getEstimator().isConvergenceAchieved()) {
			Matrix pred = new Matrix(getNumberOfObservations(), 2);
			Matrix beta = getEstimator().getParameterEstimates().getMean().getSubMatrix(0, matrixX.m_iCols - 1, 0, 0);
			Matrix linearPred = matrixX.multiply(beta);
			Matrix omega = getEstimator().getParameterEstimates().getVariance().getSubMatrix(0, matrixX.m_iCols - 1, 0, matrixX.m_iCols - 1);
			for (int i = 0; i < getNumberOfObservations(); i++) {
				pred.setValueAt(i, 0, linearPred.getValueAt(i, 0));
				Matrix x_i = matrixX.getSubMatrix(i, i, 0, matrixX.m_iCols - 1);
				pred.setValueAt(i, 1, x_i.multiply(omega).multiply(x_i.transpose()).getValueAt(0, 0));
			}
			return pred;
		} else {
			return null;
		}
	}

	@Override
	public Matrix getPredicted() {
		if (getEstimator().isConvergenceAchieved()) {
			completeLLK.setParameters(getEstimator().getParameterEstimates().getMean());
			return completeLLK.getPredictions();
		} else {
			return null;
		}
	}

	@Override
	public Matrix getResiduals() {
		if (getEstimator().isConvergenceAchieved()) {
			return y.subtract(getPredicted());
		} else {
			return null;
		}
	}

	@Override
	public CompositeLogLikelihoodWithExplanatoryVariables getCompleteLogLikelihood() {return completeLLK;}

	@Override
	public boolean isInterceptModel() {return getDataStructure().isInterceptModel();}

	@Override
	public List<String> getEffectList() {return getDataStructure().getEffectList();}

	@Override
	public double getConvergenceCriterion() {return convergenceCriterion;}
	
	public void setConvergenceCriterion(double d) {
		if (d < 1E-16) {
			throw new InvalidParameterException("The minimum value for the convergence criterion is 1E-16");
		} 
		convergenceCriterion = d;
	}



}
