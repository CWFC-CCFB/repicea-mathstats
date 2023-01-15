/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2019 Mathieu Fortin for Rouge-Epicea
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
package repicea.math.utility;

import java.security.InvalidParameterException;

/**
 * This class implements the Gamma function.
 * @author Mathieu Fortin - November 2015
 */
public class GammaUtility {
	
	private static final double[] COEF = new double[]{1.000000000000000174663, 5716.400188274341379136, -14815.30426768413909044, 14291.49277657478554025,
    		-6348.160217641458813289, 1301.608286058321874105, -108.1767053514369634679, 2.605696505611755827729, -0.7423452510201416151527e-2,
    		0.5384136432509564062961e-7, -0.4023533141268236372067e-8};		// higher precision with these coefficients

	/**
	 * This method returns the result of the Gamma function. The implementation is the Lanczos approximation. 
	 * @param z 
	 * @return a double
	 */
	public static double gamma(double z) {
		if (z <= 0d) {
			throw new InvalidParameterException("The gamma function does not support null or negative values!");
		}
	   	double result;

	    if (z < 0.5) {
	    	result = Math.PI / (Math.sin(Math.PI*z) * gamma(1-z));
	    } else {
	        z -= 1;
	        double x = COEF[0];
	        double pval;
	        for (int k = 1; k < COEF.length; k++) {
	        	pval = COEF[k];
	            x += pval/(z+k);		// i+1 = k 
	        }

	        double t = z + COEF.length - 1 - 0.5;
	        result = Math.sqrt(2d*Math.PI) * Math.pow(t, (z + .5)) * Math.exp(-t) * x;
	    }
	    return result;
	}
	
	/**
	 * This method returns the logarithm of the Gamma function.
	 * @param z
	 * @return a double
	 */
	public static double logGamma(double z) {
		return Math.log(gamma(z));
	}
	

	private static final double K = 1.461632;
	private static final double SQRT_TWICE_PI = Math.sqrt(2 * Math.PI);
	private static final double C = SQRT_TWICE_PI / Math.exp(1) - gamma(K);

	/**
	 * This method implements the approximation of the inverse Gamma function designed by 
	 * David W. Cantrell (see http://mathforum.org/kb/message.jspa?messageID=342551&tstart=0).
	 * @param z a double equal or greater than 1.
	 * @return a double
	 */
	public static double inverseGamma(double z) {
		if (z < 1) {
			throw new InvalidParameterException("The method does not accept argument smaller than 1!");
		}
		double l_x = Math.log((z + C) / SQRT_TWICE_PI);
		double w_x = lambertW(l_x/Math.exp(1));
		return l_x/w_x + .5;
	}
	
	/**
	 * This method implements the Lambert W function.
	 * Code taken from http://keithbriggs.info/software.html  
	 * @param z
	 * @return
	 */
	private static double lambertW(double z) {
		final double eps = 4e-16;
		final double em1 = 0.3678794411714423215955237701614608; 
		double p,e,t,w;
		if (z < -em1 || Double.isInfinite(z) || Double.isNaN(z)) { 
			throw new InvalidParameterException("The parameter z must be greater than -0.367879!"); 
		}
		if (z == 0d) {
			return 0d;
		}
		if (z < -em1 + 1e-4) { // series near -em1 in sqrt(q)
			double q = z + em1;
			double r = Math.sqrt(q);
			double q2=q*q;
			double q3=q2*q;
			return 	-1.0
					+2.331643981597124203363536062168*r
					-1.812187885639363490240191647568*q
					+1.936631114492359755363277457668*r*q
					-2.353551201881614516821543561516*q2
					+3.066858901050631912893148922704*r*q2
					-4.175335600258177138854984177460*q3
					+5.858023729874774148815053846119*r*q3
					-8.401032217523977370984161688514*q3*q;  // error approx 1e-16
		}
		/* initial approx for iteration... */
		if (z < 1.0) { /* series near 0 */
			p = Math.sqrt(2d * (2.7182818284590452353602874713526625*z + 1.0));
			w = -1.0 + p * (1.0 + p * (-0.333333333333333333333 + p*0.152777777777777777777777)); 
		} else {
			w = Math.log(z); /* asymptotic */
		}
		if (z > 3.0) {
			w -= Math.log(w); /* useful? */
		}
		for (int i=0; i<10; i++) { /* Halley iteration */
			e = Math.exp(w); 
			t = w*e - z;
			p = w + 1d;
			t /= e * p - 0.5 * (p + 1.0) * t/p; 
			w -= t;
			if (Math.abs(t) < eps * (1.0 + Math.abs(w))) {
				return w; /* rel-abs error */
			}
		}
		/* should never get here */
		throw new InvalidParameterException("Unable to reach convergence for z = " + z);
	}

	/**
	 * Compute the first derivative of the gamma function. <br>
	 * <br>
	 * The calculation is based on the digamma function.
	 * @see<a href=https://en.wikipedia.org/wiki/Digamma_function> Digamma function </a>
	 * @param d
	 * @return the first derivative (a double)
	 */
	public static double gammaFirstDerivative(double d) {
		return gamma(d) * digamma(d);
	}
	
	
	/**
	 * Compute an approximation of the digamma function. <br>
	 * <br>
	 * The approximation is calculated as ln(d) - 1/2d.
	 * @see<a href=https://en.wikipedia.org/wiki/Digamma_function> Digamma function </a>
	 * @param d a strictly positive double 
	 * @return a double
	 */
	public static double digamma(double d) {
		if (d <= 0d) {
			throw new InvalidParameterException("The digamma function is not defined for values smaller than or equal to 0!");
		}
		double d_star = d;
		double corrTerm = 0;
		while(d_star < 6) {
			corrTerm += 1/d_star;
			d_star += 1;
		}
		double result = getDigammaExpansion(d_star) - corrTerm;
		return result;
	}
	
	
	private static double getDigammaExpansion(double z) {
		double z2 = z*z;
		double z4 = z2*z2;
		double z6 = z4*z2;
		return Math.log(z) - 1d/(2*z) - 1d/(12*z2) + 1d/(120*z4) - 1d/(252*z4*z2) + 1d/(240*z4*z4) -
				1d/(132*z6*z4) + 691d/(32760*z6*z6) - 1d/(12*z6*z6*z2);
	}
	
	
	/**
	 * Compute an approximation of the digamma function. <br>
	 * <br>
	 * The approximation is calculated as ln(d) - 1/2d.
	 * @see<a href=https://en.wikipedia.org/wiki/Digamma_function> Digamma function </a>
	 * @param d a strictly positive double 
	 * @return a double
	 */
	public static double trigamma(double d) {
		if (d <= 0d) {
			throw new InvalidParameterException("The digamma function is not defined for values smaller than or equal to 0!");
		}
		double d_star = d;
		double corrTerm = 0;
		while(d_star < 6) {
			corrTerm += -1/(d_star*d_star);
			d_star += 1;
		}
		double expansion = getTrigammaExpansion(d_star); 
		double result = expansion - corrTerm;
		return result;
	}

	private static double getTrigammaExpansion(double z) {
		double z2 = z*z;
		double z4 = z2*z2;
		double z6 = z4*z2;
		return 1d/z + 1d/(2*z2) + 1d/(6*z2*z) - 1d/(30*z4*z) + 1d/(42*z6*z) - 1d/(30*z6*z2*z) + 
				5d/(66*z6*z4*z) - 691/(2730*z6*z6*z) + 7d/(6*z6*z6*z2*z); 
	}
	
//	public static void main(String[] args) {
//		double fake = C;
//		double z = 10;
//		double actual = GammaFunction.gamma(z);
//		double zBack = GammaFunction.inverseGamma(actual);
//		int u = 0;
//	}
	
}
