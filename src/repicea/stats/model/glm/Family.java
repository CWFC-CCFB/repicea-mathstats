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
import java.util.Arrays;
import java.util.List;

import repicea.math.MathematicalFunction;
import repicea.math.Matrix;

/**
 * The Family class defines the distribution of the response 
 * variable in a generalized linear model. 
 * 
 * @author Mathieu Fortin - October 2022
 */
public class Family {

	public static enum GLMDistribution {
		Bernoulli(0, new ArrayList<String>(), new ArrayList<Double>(), LinkFunction.Type.Logit, LinkFunction.Type.CLogLog) ,
		NegativeBinomial(1, Arrays.asList(new String[] {"Theta"}), Arrays.asList(new Double[] {1d}), LinkFunction.Type.Log);
				
		private final List<LinkFunction.Type> acceptedTypes;
		protected final int nbAdditionalParms;
		private final ArrayList<Double> startingValues;
		protected final ArrayList<String> additionalParmNames;
		
		GLMDistribution(int nbAddParms, List<String> names, List<Double> values, LinkFunction.Type... types) {
			nbAdditionalParms = nbAddParms;
			acceptedTypes = new ArrayList<LinkFunction.Type>();
			startingValues = new ArrayList<Double>();
			startingValues.addAll(values);
			additionalParmNames = new ArrayList<String>();
			additionalParmNames.addAll(names);
			if (types != null) {
				for (LinkFunction.Type t : types) {
					acceptedTypes.add(t);
				}
			}
		}
		
		boolean isAcceptedType(LinkFunction.Type t) {
			return acceptedTypes.contains(t);
		}
		
		private Matrix getAdditionalParmsStartingValues() {
			if (startingValues != null && !startingValues.isEmpty())
				return new Matrix(startingValues);
			else
				return null;
		}
		
		protected Matrix getStartingParms(int nbFixedEffectParms) {
			return getStartingParms(new Matrix(nbFixedEffectParms, 1));
		}
		
		protected Matrix getStartingParms(Matrix fixedEffectParms) {
			Matrix addParms = getAdditionalParmsStartingValues();
			return addParms != null ?
					fixedEffectParms.matrixStack(addParms, true) :
						fixedEffectParms;
		}

	}

	protected final GLMDistribution dist;
	protected final LinkFunction lf;
	
	private Family(GLMDistribution d, LinkFunction lf) {
		if (!d.isAcceptedType(lf.getType())) {
			throw new InvalidParameterException("This distribution " + d.name() + " does not accept this link function " + lf.getType().name() + "!");
		}
		dist = d;
		this.lf = lf;
	}
	
	
	
	/**
	 * Create a Family instance.
	 * 
	 * @param d a Distribution enum
	 * @param type a type of link function (see LinkFunction.Type)
	 * @param eta a MathematicalFunction (can be null)
	 * @return a Family instance
	 */
	public static Family createFamily(GLMDistribution d, LinkFunction.Type type, MathematicalFunction eta) { 
		if (d == null || type == null) {
			throw new InvalidParameterException("The arguments d and type cannot be null!");
		}
		LinkFunction lf = eta == null ? 
				new LinkFunction(type) :
					new LinkFunction(type, eta);
		return createFamily(d, lf);
	}
	
	/**
	 * Create a Family instance.
	 * 
	 * @param d a Distribution enum
	 * @param lf a link function (see LinkFunction)
	 * @return a Family instance
	 */
	public static Family createFamily(GLMDistribution d, LinkFunction lf) { 
		return new Family(d, lf);
	}

}
