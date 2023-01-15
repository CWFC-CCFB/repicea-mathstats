/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2015 Mathieu Fortin for Rouge-Epicea
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
package repicea.stats.distributions;

import java.security.InvalidParameterException;

import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.serial.xml.PostXmlUnmarshalling;
import repicea.stats.Distribution;
import repicea.stats.StatisticalUtility;
import repicea.stats.distributions.utility.GaussianUtility;

@SuppressWarnings("serial")
public class StandardGaussianDistribution implements ContinuousDistribution, PostXmlUnmarshalling {

	private static StandardGaussianDistribution Singleton;
	
	private Matrix mu;
	private Matrix sigma2;
	private Matrix lowerCholTriangle;
		
	/**
	 * This constructor creates a Gaussian distribution with mean mu 0 and variance 1.
	 */
	protected StandardGaussianDistribution() {
		Matrix mu = new Matrix(1,1);
		SymmetricMatrix sigma2 = new SymmetricMatrix(1);
		sigma2.setValueAt(0, 0, 1d);
		setMean(mu);
		setVariance(sigma2);
	}
	
	/**
	 * This method returns the single instance of the StandardGaussianDistribution class.
	 * @return a StandardGaussianDistribution instance
	 */
	public static StandardGaussianDistribution getInstance() {
		if (Singleton == null) {
			Singleton = new StandardGaussianDistribution();
		}
		return Singleton;
	}
	
	@Override
	public boolean isMultivariate() {
		return getMu().m_iRows > 1;
	}
	
	@Override
	public Matrix getRandomRealization() {
		Matrix mean = getMean();
		Matrix standardDeviation = getStandardDeviation();
		Matrix normalStandardDeviates = StatisticalUtility.drawRandomVector(standardDeviation.m_iRows, Distribution.Type.GAUSSIAN);
		return mean.add(standardDeviation.multiply(normalStandardDeviates));
	}

	/**
	 * This method returns the lower triangle of the Cholesky decomposition of the variance-covariance matrix.
	 * @return a Matrix instance
	 */
	public Matrix getStandardDeviation() {
		if (lowerCholTriangle == null) {
			lowerCholTriangle = getSigma2().getLowerCholTriangle();
		}
		return lowerCholTriangle;
	}
	
	@Override
	public Matrix getMean() {return getMu();}

	@Override
	public SymmetricMatrix getVariance() {return getSigma2();}

	@Override
	public Type getType() {return Distribution.Type.GAUSSIAN;}

	protected void setMean(Matrix mu) {
		this.mu = mu;
	}

	protected void setVariance(SymmetricMatrix sigma2) {
		if (sigma2 != null && !sigma2.isSymmetric()) {
			throw new InvalidParameterException("The variance-covariance matrix must be symmetric!");
		}
		this.sigma2 = sigma2;
		lowerCholTriangle = null;
	}

	protected Matrix getMu() {return mu;}
	
	protected SymmetricMatrix getSigma2() {return (SymmetricMatrix) sigma2;}
	
	
	@Override
	public boolean isParametric() {return true;}

	
	/**
	 * This method returns the result of the probability density function of the distribution parameter.
	 * @param yValues a single double value or a Matrix instance
	 * @return a double
	 */
	@Override
	public double getProbabilityDensity(Matrix yValues) {
		if (yValues == null || !yValues.isTheSameDimension(getMu())) {
			throw new UnsupportedOperationException("Vector y is either null or its dimensions are different from those of mu!");
		} else {
			if (!isMultivariate()) {
				double y = yValues.getValueAt(0, 0);
				double mu = getMu().getValueAt(0, 0);
				double variance = getSigma2().getValueAt(0, 0);
				return GaussianUtility.getProbabilityDensity(y, mu, variance);
			} else {
				int k = yValues.m_iRows;
				Matrix residuals = yValues.subtract(getMu());
				Matrix invSigma2 = getSigma2().getInverseMatrix();
				return 1d / (Math.pow(2 * Math.PI, 0.5 * k) * Math.sqrt(getSigma2().getDeterminant())) * Math.exp(- 0.5 * residuals.transpose().multiply(invSigma2).multiply(residuals).getSumOfElements());
			}
		}
	}

	@Override
	public void postUnmarshallingAction() {
		if (!(sigma2 instanceof SymmetricMatrix)) {
			setVariance(SymmetricMatrix.convertToSymmetricIfPossible(sigma2));	// MF20220824 to ensure a proper deserialization of previously saved meta models.
		}
	}

	
//	@Override
//	public List<double[]> getQuantile(List<double[]> probabilities) {
//		if (probabilities == null || probabilities.isEmpty()) {
//			throw new InvalidParameterException("The probabilities parameter is null or empty!");
//		} else if (isMultivariate() && probabilities.size() != ((Matrix) getMean()).m_iRows) {
//			throw new InvalidParameterException("The number of values does not correspond to the dimension of the distribution!");
//		} else {
//			List<double[]> output = new ArrayList<double[]>();
//			if (!isMultivariate()) {
//				double[] probabilityLevels = probabilities.get(0);
//				double[] quantiles = new double[probabilityLevels.length];
//				for (int i = 0; i < probabilityLevels.length; i++) {
//					double standardizedValue = (values[0] - getMean().m_afData[0][0]) / Math.sqrt(getVariance().m_afData[0][0]);
//				}
////				return GaussianUtility.getCumulativeProbability(standardizedValue);
////			} else if (values.length == 2) {
////				double std1 = Math.sqrt(getVariance().m_afData[0][0]);
////				double standardizedValue1 = (values[0] - getMean().m_afData[0][0]) / std1;
////				double std2 = Math.sqrt(getVariance().m_afData[1][1]);
////				double standardizedValue2 = (values[1] - getMean().m_afData[1][0]) / std2;
////				double correlation = getVariance().m_afData[0][1] / (std1 * std2);
////				return GaussianUtility.getBivariateCumulativeProbability(standardizedValue1, standardizedValue2, correlation);
////			}
////		}
////		return -1;
//		}
//		return null;
//	}

	
//	public static void main(String[] args) {
//		GaussianFunction gf = new GaussianFunction();
//		gf.setParameterValue(FunctionParameter.Mu, new Matrix(2,1));
//		gf.setParameterValue(FunctionParameter.Sigma2, Matrix.getIdentityMatrix(2));
//		gf.setVariableValue(FunctionVariable.yVector, new Matrix(2,1));
//		double allo = gf.getValue();
//	}

}
