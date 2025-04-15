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
package repicea.stats;

import java.io.Serializable;

import repicea.math.AbstractMatrix;
import repicea.util.REpiceaTranslator;
import repicea.util.REpiceaTranslator.TextableEnum;



/**
 * The Distribution interface provides the basic function for the GaussianDistribution, StudentTDistribution and other 
 * classes of the repicea.stat.distribution package. The interface provides the getRandomRealization() method which
 * returns a random deviate of the distribution. For the implementation of probability density (or mass) functions, see
 * the repicea.stat.distribution.utility package.
 * @author Mathieu Fortin - August 2012
 * 
 * @param <M> an AbstractMatrix-derived class that stands for the mean
 * @param <V> an AbstractMatrix-derived class that stands for the variance
 */
public interface Distribution<M extends AbstractMatrix, V extends AbstractMatrix> extends MomentGettable<M,V>, Serializable {

	public enum Type implements TextableEnum {
		GAUSSIAN("Gaussian", "Gaussienne"), 
		UNIFORM("Uniform", "Uniforme"), 
		NONPARAMETRIC("Non parametric", "Non param\u00E9trique"), 
		UNKNOWN("Unknown", "Inconnue"), 
		CHI_SQUARE("Chi squared", "Chi carr\u00E9"),
		WISHART("Wishart", "Wishart"),
		STUDENT("Student's t", "t de Student");

		
		Type(String englishText, String frenchText) {
			setText(englishText, frenchText);
		}
		
		@Override
		public void setText(String englishText, String frenchText) {
			REpiceaTranslator.setString(this, englishText, frenchText);
		}

		@Override
		public String toString() {return REpiceaTranslator.getString(this);}
	
	}
	
	/**
	 * This method returns true if the distribution is parametric or false otherwise.
	 * @return a boolean
	 */
	public boolean isParametric();
	
	
	/**
	 * This method returns true if the GaussianFunction instance is multivariate.
	 * @return a boolean
	 */
	public boolean isMultivariate();
	
	
	/**
	 * Returns true if the GaussianFunction instance is univariate.
	 * @return a boolean
	 */
	public default boolean isUnivariate() {
		return !isMultivariate();
	}
	
	/** 
	 * This method returns the type of the distribution.
	 * @return a Type enum
	 */
	public Type getType();

	
	/**
	 * This method draws a random realization from the distribution.
	 * @return the observation in a Matrix instance
	 */
	public M getRandomRealization();


}
