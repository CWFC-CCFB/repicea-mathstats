/*
 * This file is part of the repicea-mathstats library.
 *
 * Copyright (C) 2024 His Majesty the King in Right of Canada
 * Author: Mathieu Fortin, Canadian Wood Fibre Centre, Canadian Forest Service
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

import java.security.InvalidParameterException;

import repicea.math.ComplexMatrix;
import repicea.math.ComplexNumber;
import repicea.math.ComplexSymmetricMatrix;
import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.math.utility.ErrorFunctionUtility;
import repicea.math.utility.GaussianUtility;
import repicea.stats.Distribution;
import repicea.stats.StatisticalUtility;
import repicea.stats.estimates.ComplexMonteCarloEstimate;
import repicea.stats.estimates.MonteCarloEstimate;
import repicea.stats.estimators.OLSEstimator;

/**
 * A class implementing various estimators of the logarithmic back transformation of
 * linear model predictions.
 * @author Mathieu Fortin - September 2024
 */
public class LogBackTransformation {

	public enum Estimator {
		/**
		 * In the context of a pure linear model, this estimator is that of Baskerville.
		 */
		Naive,
		BeauchampAndOlson,
		MonteCarlo,
		ComplexMonteCarlo
	}
	
	private static int  NB_INNER_REALIZATIONS = 10000;
	static double VERY_SMALL = 1E-8;


	/**
	 * Provide the mean predicted values on the original scale. <p>
	 * The method assumes that the response variable has been transformed using 
	 * a log transformation: w = ln(y + c), where c is a transformation offset, 
	 * typically 1. <p>
	 * 
	 * The method returns a column vector that contains the mean predicted values.
	 * If the method argument is set to Estimator.ComplexMonteCarlo, the method
	 * returns a two-column matrix, with the second column containing the variance
	 * of the mean predicted values. 
	 * 
	 * @param model a LinearModel instance
	 * @param transformationOffset the transformation offset
	 * @param method an Enum estimator - the estimator to be used
	 * @return a Matrix instance
	 */
	public static Matrix getMeanPredictedValuesOnOriginalScale(LinearModel model, double transformationOffset, Estimator method) {
		return getMeanPredictedValuesOnOriginalScale(model, null, transformationOffset, method);
	}

	
	/**
	 * Provide the mean predicted values on the original scale. <p>
	 * The method assumes that the response variable has been transformed using 
	 * a log transformation: w = ln(y + c), where c is a transformation offset, 
	 * typically 1. <p>
	 * 
	 * The method returns a column vector that contains the mean predicted values.
	 * If the method argument is set to Estimator.ComplexMonteCarlo, the method
	 * returns a two-column matrix, with the second column containing the variance
	 * of the mean predicted values. 
	 * 
	 * @param model a LinearModel instance
	 * @param xMat a design matrix for which the predictions are needed
	 * @param transformationOffset the transformation offset
	 * @param method an Enum estimator - the estimator to be used
	 * @return a Matrix instance
	 */
	public static Matrix getMeanPredictedValuesOnOriginalScale(LinearModel model, Matrix xMat, double transformationOffset, Estimator method) {
		if (model instanceof LinearModelWithTruncatedGaussianErrorTerm) {
			double sigma2 = model.getResidualVariance();
			double sigma = Math.sqrt(sigma2);
			Matrix xBeta = ((LinearModelWithTruncatedGaussianErrorTerm) model).getXBeta(xMat);
			Matrix meanValues = new Matrix(xBeta.m_iRows, 1);
			double truncation = ((LinearModelWithTruncatedGaussianErrorTerm) model).truncation;
			for (int i = 0; i < xBeta.m_iRows; i++) {
				double xBeta_i = xBeta.getValueAt(i, 0);
				double F_t = GaussianUtility.getCumulativeProbability((truncation - xBeta_i)/sigma);
				double meanValue =  (F_t < VERY_SMALL) ?
						Math.exp(xBeta_i + 0.5 * sigma2) :
							Math.exp(xBeta_i + 0.5 * sigma2) * ((1 + ErrorFunctionUtility.erf((xBeta_i + sigma2 - truncation)/Math.sqrt(2*sigma2)))/(2*(1-F_t))); 
				meanValues.setValueAt(i, 0, meanValue - transformationOffset);
			}
			return meanValues;
		} 
		
//		// Otherwise it is a simple LinearModel instance
//		if (!(model.getEstimator() instanceof OLSEstimator)) {
//			throw new UnsupportedOperationException("This method only works with models fitted using an OLS estimator!");
//		}
		Matrix pred;
		Matrix omegaChol;
		
		double sigma2hat = model.getResidualVariance();
		
		switch(method) {
		case Naive:
			pred = model.getPredicted(xMat).scalarAdd(0.5 * sigma2hat).expMatrix().scalarAdd(-transformationOffset);
			return pred;
		case BeauchampAndOlson:
			Matrix xBeta = model.getPredicted(xMat);
			pred = new Matrix(xMat.m_iRows, 1);
			SymmetricMatrix invXtX = ((OLSEstimator) model.getEstimator()).getInverseXtXMatrix();
			int n = model.getNumberOfObservations();
			for (int ii = 0; ii < xMat.m_iRows; ii++) {
				double term1 = Math.exp(xBeta.getValueAt(ii, 0) + sigma2hat * .5);
				Matrix xMat_ii = xMat.getSubMatrix(ii, ii, 0, xMat.m_iCols - 1);
				Matrix xInvXtXxT = xMat_ii.multiply(invXtX).multiply(xMat_ii.transpose());
				double phi = xInvXtXxT.getValueAt(0, 0) * n;
				double sigma2hat2 = sigma2hat * sigma2hat;
				double term2 = 1 - sigma2hat * (2*phi + sigma2hat)/(4*n) + 
						sigma2hat2 * (sigma2hat2 + 2 * (16d/3 + 2 * phi) * sigma2hat + 4 * phi * phi + 16 * phi) / (32 * n * n);
				pred.setValueAt(ii, 0, term1 * term2 - transformationOffset);
			}
			return pred;
		case MonteCarlo:
			MonteCarloEstimate mcEstimator = new MonteCarloEstimate();
			omegaChol = ((OLSEstimator) model.getEstimator()).getInverseXtXMatrix().getLowerCholTriangle();
			for (int innerReal = 0; innerReal < NB_INNER_REALIZATIONS; innerReal++) {
				mcEstimator.addRealization(getRandomDeviateFromModel(model, xMat, omegaChol)); 
			}
			return mcEstimator.getMean().scalarAdd(-transformationOffset);
		case ComplexMonteCarlo:
			ComplexMonteCarloEstimate cmcEstimator = new ComplexMonteCarloEstimate();
			omegaChol = ((OLSEstimator) model.getEstimator()).getInverseXtXMatrix().getLowerCholTriangle();
			for (int innerReal = 0; innerReal < NB_INNER_REALIZATIONS; innerReal++) {
				ComplexMatrix complexRealizations = getComplexRandomDeviateFromModel(model, xMat, omegaChol);
				cmcEstimator.addRealization(complexRealizations); 
			}
			ComplexMatrix mean = cmcEstimator.getMean();
			ComplexSymmetricMatrix cmcPsVar = cmcEstimator.getPseudoVariance();
			pred = new Matrix(mean.m_iRows, 2);
			for (int i = 0; i < mean.m_iRows; i++) {
				pred.setValueAt(i, 0, mean.getValueAt(i, 0).realPart - transformationOffset);
				pred.setValueAt(i, 1, -cmcPsVar.getValueAt(i, i).realPart);
			}
			return pred;
		default:
			throw new InvalidParameterException("The following estimator is not supported yet: " + method.name());
		}
		
	}
	
	private static Matrix getRandomDeviateFromModel(LinearModel model, Matrix xMat, Matrix omegaChol) {
		int upsilon = ((OLSEstimator) model.getEstimator()).getResidualVariance().getDegreesOfFreedom();
		double chiSquareDeviate = StatisticalUtility.getRandom().nextChiSquare(upsilon);
		double sigma2hat = model.getResidualVariance();
		double sigma2_b = sigma2hat + sigma2hat * (chiSquareDeviate/upsilon - 1);
		Matrix betaHat = model.getParameters();
		Matrix betaHat_b = betaHat.add(omegaChol.multiply(StatisticalUtility.drawRandomVector(betaHat.m_iRows, Distribution.Type.GAUSSIAN)).scalarMultiply(Math.sqrt(sigma2_b)));
		Matrix result = xMat.multiply(betaHat_b).scalarAdd(0.5 * sigma2_b);
		return result.expMatrix();
	}
	
	private static ComplexMatrix getComplexRandomDeviateFromModel(LinearModel model, Matrix xMat, Matrix omegaChol) {
		int upsilon = ((OLSEstimator) model.getEstimator()).getResidualVariance().getDegreesOfFreedom();
		double chiSquareDeviate = StatisticalUtility.getRandom().nextChiSquare(upsilon);
		double sigma2hat = model.getResidualVariance();
		ComplexNumber sigma2Hat_b = new ComplexNumber(sigma2hat, sigma2hat * (chiSquareDeviate/upsilon - 1));
		ComplexNumber sigmaHat_b = sigma2Hat_b.sqrt();
		Matrix betaHat = model.getParameters();
		Matrix betaDeviatesMat = omegaChol.multiply(StatisticalUtility.drawRandomVector(betaHat.m_iRows, Distribution.Type.GAUSSIAN)); // here we get C * epsilon
		
		ComplexNumber[] cnArrayMean = new ComplexNumber[xMat.m_iRows];
		for (int ii = 0; ii < xMat.m_iRows; ii++) {
			Matrix xMat_ii = xMat.getSubMatrix(ii, ii, 0, xMat.m_iCols - 1);
			ComplexNumber betaDeviates = new ComplexNumber(0, xMat_ii.multiply(betaDeviatesMat).getValueAt(0, 0)); // here we get 0 + x C epsilon i 
			betaDeviates = betaDeviates.multiply(sigmaHat_b); // here we get 0 + sqrt{sigma2Hat_b} x C epsilon i
			double xBetaHat = xMat_ii.multiply(betaHat).getValueAt(0, 0); 
			ComplexNumber mean = betaDeviates.add(xBetaHat); // here we get x beta + + sqrt{sigma2Hat_b} x C epsilon i
			ComplexNumber meanPlusCorrectionFactor = mean.add(sigma2Hat_b.multiply(0.5)); // here we get x beta + + sqrt{sigma2Hat_b} x C epsilon i + 0.5 * sigma2_b  				
			cnArrayMean[ii] = meanPlusCorrectionFactor.exp();
		}
		return new ComplexMatrix(cnArrayMean); 
	}

	/**
	 * Set the number of realizations to be performed for Monte Carlo estimators.
	 * @param nbRealizations an integer ranging between 1 and 1E7
	 */
	public static void setInnerRealizationsForMonteCarloEstimators(int nbRealizations) {
		if (nbRealizations < 1 || nbRealizations > 1E7) {
			throw new InvalidParameterException("The nbRealizations argument must be an integer ranging between 1 and 1E7 !");
		}
		NB_INNER_REALIZATIONS = nbRealizations;
	}
	

	/**
	 * Provide the residual variances on the original scale. <p>
	 * The method assumes that the response variable has been transformed using 
	 * a log transformation: w = ln(y + c), where c is a transformation offset, 
	 * typically 1. <p>
	 * 
	 * The method returns a column vector that contains the residual variances. No 
	 * offset is required since it only affects the mean and not the variance. In its
	 * current form, the method only support the Estimator.Naive enum.
	 * 
	 * @param model a LinearModel instance
	 * @param method an Enum estimator - the estimator to be used
	 * @return a Matrix instance
	 */
	public static Matrix getResidualVariancesOnOriginalScale(LinearModel model, Estimator method) {
		return 	getResidualVariancesOnOriginalScale(model, null, method);
	}

	
	
	/**
	 * Provide the residual variances on the original scale. <p>
	 * The method assumes that the response variable has been transformed using 
	 * a log transformation: w = ln(y + c), where c is a transformation offset, 
	 * typically 1. <p>
	 * 
	 * The method returns a column vector that contains the residual variances. No 
	 * offset is required since it only affects the mean and not the variance. In its
	 * current form, the method only support the Estimator.Naive enum.
	 * 
	 * @param model a LinearModel instance
	 * @param xMat a design matrix for which the predictions are needed
	 * @param method an Enum estimator - the estimator to be used
	 * @return a Matrix instance
	 */
	public static Matrix getResidualVariancesOnOriginalScale(LinearModel model, Matrix xMat, Estimator method) {
		if (model instanceof LinearModelWithTruncatedGaussianErrorTerm) {
			double sigma2 = model.getResidualVariance();
			double sigma = Math.sqrt(sigma2);
			Matrix xBeta = ((LinearModelWithTruncatedGaussianErrorTerm) model).getXBeta(xMat);
			Matrix residualVariances = new Matrix(xBeta.m_iRows, 1);
			double truncation = ((LinearModelWithTruncatedGaussianErrorTerm) model).truncation;
			for (int i = 0; i < xBeta.m_iRows; i++) {
				double xBeta_i = xBeta.getValueAt(i, 0);
				double F_t = GaussianUtility.getCumulativeProbability((truncation - xBeta_i)/sigma);
				double meanValue =  (F_t < VERY_SMALL) ?
						Math.exp(xBeta_i + 0.5 * sigma2) :
							Math.exp(xBeta_i + 0.5 * sigma2) * ((1 + ErrorFunctionUtility.erf((xBeta_i + sigma2 - truncation)/Math.sqrt(2*sigma2)))/(2*(1-F_t))); 
				double sqrt2Sigma2 = Math.sqrt(2 * sigma2);
				double part1 = Math.exp(2*xBeta_i + 2* sigma2) * (1 + ErrorFunctionUtility.erf((xBeta_i + 2*sigma2 - truncation)/sqrt2Sigma2));
				double part2 = -2 * meanValue * Math.exp(xBeta_i + .5*sigma2) * (1 + ErrorFunctionUtility.erf((xBeta_i + sigma2 - truncation)/sqrt2Sigma2));
				double part3 = meanValue * meanValue * (1 - ErrorFunctionUtility.erf((truncation - xBeta_i)/sqrt2Sigma2));
				double var = 1d / (2 - 2*F_t) * (part1 + part2 + part3);
				residualVariances.setValueAt(i, 0, var);
			}
			return residualVariances;
		} 
		
//		// Otherwise it is a simple LinearModel instance
//		if (!(model.getEstimator() instanceof OLSEstimator)) {
//			throw new UnsupportedOperationException("This method only works with models fitted using an OLS estimator!");
//		}
		
		double sigma2hat = model.getResidualVariance();
		
		switch(method) {
		case Naive:
			Matrix var = model.getPredicted(xMat).scalarMultiply(2d).scalarAdd(sigma2hat).expMatrix().scalarMultiply(Math.exp(sigma2hat) - 1);
			return var;
		default:
			throw new InvalidParameterException("The following estimator is not supported yet: " + method.name());
		}
		
	}
	

	
}
