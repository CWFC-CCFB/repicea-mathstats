/*
 * This file is part of the repicea-statistics library.
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
package repicea.stats.model.glm.copula;

import java.security.InvalidParameterException;

import repicea.math.Matrix;
import repicea.math.optimizer.AbstractOptimizer.OptimizationException;
import repicea.stats.data.DataSet;
import repicea.stats.data.GenericHierarchicalSpatialDataStructure;
import repicea.stats.data.HierarchicalStatisticalDataStructure;
import repicea.stats.data.StatisticalDataException;
import repicea.stats.model.CompositeLogLikelihoodWithExplanatoryVariables;
import repicea.stats.model.glm.Family.GLMDistribution;
import repicea.stats.model.glm.GeneralizedLinearModel;

/**
 * The FGMCopulaGLModel class implements a Farlie-Gumbel-Morgenstein (FGM) copula for binary outcomes. The method is described in 
 * Bhat and Sener 2009.
 * @see <a href=http://www.springerlink.com/content/t33474476r231w60/> 
 * Bhat, C.R., and Sener, I. N. 2009. A copula-based closed form binary logit choice model for accomodating spatial correlation across observational units. Journal of Geographical Systems 11(3): 243-272.
 * </a>
 * @author Mathieu Fortin - June 2011
 */
public class FGMCopulaGLModel extends GeneralizedLinearModel {
	
	private final CopulaExpression copula;
	
	/**
	 * Constructor.
	 * @param glm a GeneralizedLinearModel instance
	 * @param copula a CopulaExpression instance
	 * @throws StatisticalDataException if the hierarchical level specification in the copula is not found in the data set
	 * @throws OptimizationException if the preliminary generalized linear model fails to converge
	 */
	public FGMCopulaGLModel(GeneralizedLinearModel glm, CopulaExpression copula) throws StatisticalDataException, OptimizationException {
		super(glm);
		if (glm.getDistribution() != GLMDistribution.Bernoulli) {
			throw new InvalidParameterException("The FGM copula only supports the Bernoulli distribution!");
		}
		if (!glm.getEstimator().isConvergenceAchieved()) {
			glm.doEstimation();
		} else {
			throw new OptimizationException("The generalized linear model could not be fitted!");
		}
		glm.setParameters(glm.getEstimator().getParameterEstimates().getMean());
		this.copula = copula;
		this.copula.initialize(this, getDataStructure());
		getCompleteLogLikelihood().initialize(getDataStructure(), copula);
	}
	
	@Override
	public FGMCompositeLogLikelihood getCompleteLogLikelihood() {
		return (FGMCompositeLogLikelihood) super.getCompleteLogLikelihood();
	}
	
	@Override
	public Matrix getParameters() {return individualLLK.getParameters().matrixStack(copula.getParameters(), true);}
	

	@Override
	public void setParameters(Matrix beta) {
		if (copula == null) {
			individualLLK.setParameters(beta);
		} else {
			individualLLK.setParameters(beta.getSubMatrix(0, beta.m_iRows - copula.getNumberOfParameters() - 1, 0, 0));
			copula.setParameters(beta.getSubMatrix(beta.m_iRows - copula.getNumberOfParameters(), beta.m_iRows - 1, 0, 0));
		}
//		if (beta == null) {
//			individualLLK.setParameters(new Matrix(matrixX.m_iCols, 1));		// default starting parameters at 0
//		} else {
//			individualLLK.setParameters();
//			copula.setParameters(beta.getSubMatrix(beta.m_iRows - copula.getNumberOfParameters(), beta.m_iRows - 1, 0, 0));
//		}
	}

	@Override
	protected CompositeLogLikelihoodWithExplanatoryVariables createCompleteLLK(Object addParm) {
		return new FGMCompositeLogLikelihood(individualLLK,	matrixX, y);
	}
	
	@Override
	public String toString() {
		return "Generalized linear model based on FGM copula";
	}

//	@Override
//	public String getSummary() {
//		return super.getSummary();
//		if (getEstimator().isConvergenceAchieved()) {
//			NumberFormat formatter = NumberFormat.getInstance();
//			formatter.setMaximumFractionDigits(4);
////			SpearmanCorrelationCoefficient scc = new SpearmanCorrelationCoefficient();
////			CorrelationEstimate[] spearmanCorr = scc.getSpearmanCorrelationCoefficient(this);
////			System.out.println("Spearman's correlation coefficients");
////			for (int i = 0; i < spearmanCorr.length; i++) {
////				CorrelationEstimate estimate = spearmanCorr[i];
////				System.out.println(i + " : " + "Estimate : " + formatter.format(estimate.getMean()) +
////						"; n : " + estimate.getSampleSize() + 
////						"; t-value : " + formatter.format(estimate.getStudentT()));
////			}
//		}
//	}
	
	protected CopulaExpression getCopula() {return copula;}
	
	
	@Override
	protected HierarchicalStatisticalDataStructure createDataStructure(DataSet dataSet, Object addParm) {
		return new GenericHierarchicalSpatialDataStructure(dataSet);
	}

	@Override
	protected HierarchicalStatisticalDataStructure getDataStructure() {
		return (HierarchicalStatisticalDataStructure) super.getDataStructure();
	}

	
}
