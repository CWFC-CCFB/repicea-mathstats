/*
 * This file is part of the repicea-mathstats library.
 *
 * Copyright (C) 2023 His Majesty the King in Right of Canada
 * Author: Mathieu Fortin, Canadian Wood Fibre Centre, CFS
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
package repicea.stats.model.lm;

import repicea.math.AbstractMathematicalFunction;
import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.math.utility.ErrorFunctionUtility;
import repicea.math.utility.GaussianUtility;
import repicea.stats.data.DataSet;
import repicea.stats.model.CompositeLogLikelihood;
import repicea.stats.model.CompositeLogLikelihoodWithExplanatoryVariables;

/**
 * A class that implements the truncated TOBIT model.
 * @see <a href=http://doi.org/https://doi.org/10.1093/forestry/cpae017> 
 * Fortin, M. 2024. The impact of natural constraints in linear regression of log transformed response variables. Forestry, cpae017</a>
 */
@SuppressWarnings("serial")
public class LinearModelWithTruncatedGaussianErrorTerm extends LinearModel {

	class TruncatedGaussianIndividualLogLikelihood extends GaussianIndividualLogLikelihood {

		class TruncationFunction extends AbstractMathematicalFunction {

			double getCumulativeProbability() {
				double xBeta = TruncatedGaussianIndividualLogLikelihood.this.getOriginalFunction().getValue();
				double standardizedThreshold = (truncation - xBeta) / Math.sqrt(sigma2);
				return GaussianUtility.getCumulativeProbability(standardizedThreshold);
			}
			
			@Override
			public Double getValue() {
				return - Math.log(1 - getCumulativeProbability());
			}

			private Matrix get_dF_dAll() {
				double xBeta = TruncatedGaussianIndividualLogLikelihood.this.getOriginalFunction().getValue();
				double f_t = GaussianUtility.getProbabilityDensity(truncation, xBeta, sigma2);
				Matrix xVect = TruncatedGaussianIndividualLogLikelihood.this.getOriginalFunction().getGradient();
				Matrix dF_dB = xVect.scalarMultiply( -1 * f_t);
				double res = truncation - xBeta;
				double F_t = getCumulativeProbability();
				double erfParameter = res / Math.sqrt(2 * sigma2);
				double dF_ds = -1/(2*sigma2) * F_t + (ErrorFunctionUtility.erf(erfParameter) + 1)/(4*sigma2) - f_t * res / (2*sigma2);
				Matrix dF_dAll = new Matrix(TruncatedGaussianIndividualLogLikelihood.this.getNumberOfParameters(),1);
				dF_dAll.setSubMatrix(dF_dB, 0, 0);
				dF_dAll.setValueAt(dF_dAll.m_iRows - 1, 0, dF_ds);
				return dF_dAll;
			}
			
			@Override
			public Matrix getGradient() {
				Matrix dF_dAll = get_dF_dAll();
				return dF_dAll.scalarMultiply(1 / (1 - getCumulativeProbability()));
			}

			@Override
			public SymmetricMatrix getHessian() {
				Matrix dF_dAll = get_dF_dAll();
				double xBeta = TruncatedGaussianIndividualLogLikelihood.this.getOriginalFunction().getValue();
				double res = truncation - xBeta;
				double erfParameter = res / Math.sqrt(2 * sigma2);
				double f_t = GaussianUtility.getProbabilityDensity(truncation, xBeta, sigma2);
				Matrix xVect = TruncatedGaussianIndividualLogLikelihood.this.getOriginalFunction().getGradient();
				Matrix dF_dB = xVect.scalarMultiply( -1 * f_t);
				double F_t = getCumulativeProbability();
				
				int nbParms = TruncatedGaussianIndividualLogLikelihood.this.getNumberOfParameters();
				
				Matrix d2F_d2B = dF_dB.scalarMultiply((truncation - xBeta) / sigma2).multiply(xVect.transpose());
				Matrix d2F_dBds = dF_dB.scalarMultiply(-1/(2*sigma2) + (res * res / (2*sigma2*sigma2)));
				double d2F_d2s = 1 / (4*sigma2*sigma2) * 
						(3 *F_t - 1.5 * (ErrorFunctionUtility.erf(erfParameter) + 1 ) - f_t * res * (res*res/sigma2 - 3));
				double oneMinusF_t = 1 - F_t;
				Matrix d2Fd2All = new Matrix(nbParms, nbParms);
				d2Fd2All.setSubMatrix(d2F_d2B, 0, 0);
				d2Fd2All.setSubMatrix(d2F_dBds, 0, d2Fd2All.m_iCols-1);
				d2Fd2All.setSubMatrix(d2F_dBds.transpose(), d2Fd2All.m_iRows - 1, 0);
				d2Fd2All.setValueAt(d2Fd2All.m_iRows - 1, d2Fd2All.m_iRows - 1, d2F_d2s);
				
				Matrix hessianPart1 = dF_dAll.multiply(dF_dAll.transpose()).scalarMultiply(1/(oneMinusF_t * oneMinusF_t));
				Matrix hessianPart2 = d2Fd2All.scalarMultiply(1 / oneMinusF_t);
				return SymmetricMatrix.convertToSymmetricIfPossible(hessianPart1.add(hessianPart2));
			}
			
		}

		TruncationFunction tf;
		
		TruncatedGaussianIndividualLogLikelihood(int nbParmsFromMatrixX, Matrix startingValues) {
			super(nbParmsFromMatrixX, startingValues);
			tf = new TruncationFunction();
		}
		
		@Override
		public Double getValue() {
			double additionalTerm = tf.getValue();
			return super.getValue() + additionalTerm;
		}

		
		@Override
		public Matrix getGradient() {
			Matrix additionalTerm = tf.getGradient();
			return super.getGradient().add(additionalTerm);
		}

		@Override
		public SymmetricMatrix getHessian() {
			Matrix additionalTerm = tf.getHessian();
			return SymmetricMatrix.convertToSymmetricIfPossible(super.getHessian().add(additionalTerm));
		}
		
	}

	protected final double truncation;
	
	/**
	 * Constructor for maximum likelihood estimation.
	 * @param dataSet a DataSet instance
	 * @param modelDefinition a model definition
	 * @param startingValues a Matrix of starting values
	 * @param truncation a threshold (values below this threshold are inconsistent)
	 */
	public LinearModelWithTruncatedGaussianErrorTerm(DataSet dataSet, String modelDefinition, Matrix startingValues, double truncation) {
		super(dataSet, modelDefinition, startingValues);
		this.truncation = truncation;
	}
	
	@Override
	public Matrix getPredicted(Matrix xMatrix) {		
		Matrix xBeta = getXBeta(xMatrix);
		Matrix parms = getParameters();
		double sigma2 = parms.getValueAt(parms.m_iRows - 1, 0);
		double sigma = Math.sqrt(sigma2);
		Matrix pred = new Matrix(xBeta.m_iRows, 1);
		for (int i = 0; i < xBeta.m_iRows; i++) {
			double mu = xBeta.getValueAt(i, 0);
			double prediction = mu + sigma2 * GaussianUtility.getProbabilityDensity(truncation, mu, sigma2) / (1 - GaussianUtility.getCumulativeProbability((truncation - mu)/sigma));
			pred.setValueAt(i, 0, prediction);
		}
		return pred;
	}

	Matrix getXBeta(Matrix xMatrix) {
		Matrix xMat = xMatrix != null ?
				xMatrix :
					getMatrixX();
		Matrix parms = getParameters();
		Matrix beta = parms.getSubMatrix(0, xMat.m_iCols - 1, 0, 0);
		Matrix xBeta = xMat.multiply(beta);
		return xBeta;
	}
	
//	public Matrix getPredOnLogBackTransformedScale(Matrix xMatrix, double transformationOffset, boolean varianceRequired) {
//		double sigma2 = getResidualVariance();
//		double sigma = Math.sqrt(sigma2);
//		Matrix xBeta = getXBeta(xMatrix);
//		Matrix meanValues = new Matrix(xBeta.m_iRows, varianceRequired ? 2 : 1);
//		for (int i = 0; i < xBeta.m_iRows; i++) {
//			double xBeta_i = xBeta.getValueAt(i, 0);
//			double F_t = GaussianUtility.getCumulativeProbability((truncation - xBeta_i)/sigma);
//			double meanValue =  (F_t < LogBackTransformation.VERY_SMALL) ?
//					Math.exp(xBeta_i + 0.5 * sigma2) :
//						Math.exp(xBeta_i + 0.5 * sigma2) * ((1 + ErrorFunctionUtility.erf((xBeta_i + sigma2 - truncation)/Math.sqrt(2*sigma2)))/(2*(1-F_t))); 
//			meanValues.setValueAt(i, 0, meanValue - transformationOffset);
//			if (varianceRequired) {
//				double sqrt2Sigma2 = Math.sqrt(2 * sigma2);
//				double part1 = Math.exp(2*xBeta_i + 2* sigma2) * (1 + ErrorFunctionUtility.erf((xBeta_i + 2*sigma2 - truncation)/sqrt2Sigma2));
//				double part2 = -2 * meanValue * Math.exp(xBeta_i + .5*sigma2) * (1 + ErrorFunctionUtility.erf((xBeta_i + sigma2 - truncation)/sqrt2Sigma2));
//				double part3 = meanValue * meanValue * (1 - ErrorFunctionUtility.erf((truncation - xBeta_i)/sqrt2Sigma2));
//				double var = 1d / (2 - 2*F_t) * (part1 + part2 + part3);
//				meanValues.setValueAt(i, 1, var);
//			}
//		}
//		return meanValues;
//	}
	
	@Override
	public String toString() {
		return "Linear model with residual errors following a truncated Gaussian distribution (Truncated Tobit model)" + System.lineSeparator() +"fitted with maximum likelihood estimator";
	}

	
	@Override
	public CompositeLogLikelihood getCompleteLogLikelihood() {
		return new CompositeLogLikelihoodWithExplanatoryVariables(new TruncatedGaussianIndividualLogLikelihood(getMatrixX().m_iCols, startingValues), getMatrixX(), getVectorY());
	}

}
