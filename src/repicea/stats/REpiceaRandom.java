/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2017 Mathieu Fortin for Rouge-Epicea
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
package repicea.stats;

import java.security.InvalidParameterException;
import java.util.Random;

import repicea.math.Matrix;
import repicea.stats.distributions.utility.NegativeBinomialUtility;
import repicea.stats.distributions.utility.TruncatedGaussianUtility;
import repicea.stats.distributions.utility.WeibullUtility;

@SuppressWarnings("serial")
public class REpiceaRandom extends Random {
	
	private static final double OneThird = 1d/3;
	
	protected REpiceaRandom() {
		super();
	}
	
	private double getRandomGammaForShapeGreaterThanOrEqualToOne(double shape) {
		double d = shape - OneThird;
		double c = 1d / Math.sqrt(9 * d);
		boolean found = false;
		double z, u;
		double v = 0d;
		while (!found) {
			z = nextGaussian();
			u = nextDouble();
			v = Math.pow(1 + c * z, 3d);
			boolean firstCondition = z > -1d/c;
			boolean secondCondition = Math.log(u) < .5 * z*z + d - d * v + d * Math.log(v);
			if (firstCondition) {
				if (secondCondition) {
					found = true; 
				}
			}
		}
		return d * v;
	}

	private double getRandomGammaForAnyShape(double shape) {
		if (shape >= 1) {
			return getRandomGammaForShapeGreaterThanOrEqualToOne(shape);
		} else {
			double x = this.getRandomGammaForShapeGreaterThanOrEqualToOne(shape + 1);
			return x * Math.pow(nextDouble(), 1d/shape);
		}
	}

	
	/**
	 * Returns a random deviate from a beta distribution.
	 *  with scales scale1 and beta.
	 * @param scale1 a double larger than 0
	 * @param scale2 a double larger than 0
	 * @return a double
	 */
	public double nextBeta(double scale1, double scale2) {
		double x = nextGamma(scale1, 1d);
		double y = nextGamma(scale2, 1d);
		return x / (x + y);
	}
	
	
	/**
	 * Returns a random gamma distributed value following Marsaglia and Tsang's method. The 
	 * mean of the distribution is obtained through the product of the shape and the scale.
	 * @param shape a double larger than 0
	 * @param scale a double larger than 0
	 * @return a double
	 */
	public double nextGamma(double shape, double scale) {
		if (shape <= 0d || scale <= 0d) {
			throw new InvalidParameterException("The shape and the scale must be larger than 0!");
		}
		double x = getRandomGammaForAnyShape(shape);
		return x * scale;
	}
	
	
	/**
	 * This method returns a random integer that follows negative binomial distribution.
	 * @param mean the mean of the distribution
	 * @param dispersion the dispersion parameter
	 * @return an integer
	 */
	public int nextNegativeBinomial(double mean, double dispersion) {
		double threshold = nextDouble();	// to determine how many recruits there are
		double cumulativeProb = 0.0;
		int output = -1;
		
		while (threshold > cumulativeProb) {		
			output++;
			double massProb = NegativeBinomialUtility.getMassProbability(output, mean, dispersion);
			cumulativeProb += massProb;
		}
		return output;
	}

	
	/**
	 * Returns a random deviate from the standard Student's t distribution. The algorithm behind 
	 * the random deviate generation is that of Bailey (1994) based on polar generation.
	 * @param df the degrees of freedom
	 * @see Bailey, R.W. 1994. Polar generation of random variances with the t-distribution. 
	 * Mathematics of Computation 62(206): 779-781.
	 */
	public double nextStudentT(double df) {
		double W = 2d;
		double U = 0;
		while (W > 1) {
			U = nextDouble();
			double V = nextDouble();
			U = 2 * U - 1;
			V = 2 * V - 1;
			W = U * U + V * V;
		}
		double C2 = U * U / W;
		double R2 = df * (Math.pow(W, - 2d / df) - 1);
		double result;
		if (nextDouble() < .5) {
			result = Math.sqrt(R2*C2);
		} else {
			result = - Math.sqrt(R2*C2);
		}
		return result;
	}

	
	/**
	 * This method returns a Chi squared random value.
	 * @param df the degrees of freedom
	 * @return a double
	 */
	public double nextChiSquare(int df) {
		if (df <= 0) {
			throw new InvalidParameterException("The number of degrees of freedom should be larger than 0");
		}
		double sumSquared = 0;
		for (int i = 0; i < df; i++) {
			double gaussian = nextGaussian();
			sumSquared += gaussian * gaussian;
		}
		return sumSquared;
	}

	/**
	 * This method returns the matrix A in the Bartlett decomposition.
	 * @param df degrees of freedom
	 * @param dim the dimensions of the matrix
	 * @return a Matrix
	 */
	public Matrix nextBartlettDecompositionMatrix(int df, int dim) {
		Matrix aMat = new Matrix(dim, dim);
		for (int i = 0; i < aMat.m_iRows; i++) {
			for (int j = 0; j <= i; j++) {
				if (i == j) {
					aMat.setValueAt(i, j, Math.sqrt(nextChiSquare(df - i)));	
				} else {
					aMat.setValueAt(i, j, nextGaussian());
				}
			}
		}
		return aMat;
	}

	/**
	 * Return a random realization of a Weibull distribution. <br>
	 * <br>
	 * The method first draws a realization from a uniform distribution [0,1] and
	 * then computes the quantile.
	 * 
	 * @param k the shape parameter (must be greater than 0)
	 * @param lambda the scale parameter (must be greater than 0)
	 * @param theta the location parameter
	 * @return the realization
	 */
	public double nextWeibull(double k, double lambda, double theta) {
		return WeibullUtility.getQuantile(nextDouble(), k, lambda, theta);
	}
	
	/**
	 * Return a random realization of a Weibull distribution without location parameter. <br>
	 * <br>
	 * The method first draws a realization from a uniform distribution [0,1] and
	 * then computes the quantile.
	 * 
	 * @param k the shape parameter (must be greater than 0)
	 * @param lambda the scale parameter (must be greater than 0)
	 * @return the realization
	 */
	public double nextWeibull(double k, double lambda) {
		return nextWeibull(k, lambda, 0d);
	}

	/**
	 * Return a random realization of a truncated Gaussian distribution. <br>
	 * <br>
	 * The method first draws a realization from a uniform distribution [0,1] and
	 * then computes the quantile.
	 * 
	 * @param mu the mean of the original normal distribution
	 * @param sigma2 the variance of the original distribution
	 * @param a the lower bound
	 * @param b the upper bound
	 * @return the realization
	 */
	public double nextTruncatedGaussian(double mu, double sigma2, double a, double b) {
		return TruncatedGaussianUtility.getQuantile(nextDouble(), mu, sigma2, a, b);
	}
}
