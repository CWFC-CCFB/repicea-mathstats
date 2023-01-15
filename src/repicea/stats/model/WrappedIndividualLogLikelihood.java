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
package repicea.stats.model;

import repicea.math.LogFunctionWrapper;
import repicea.math.Matrix;


/**
 * The WrappedIndividualLogLikelihood wrapper a likelihood function in a log wrapper in order 
 * to compute the log-likelihood of individual observations.
 * @author Mathieu Fortin - June 2022
 */
@SuppressWarnings("serial")
public class WrappedIndividualLogLikelihood extends LogFunctionWrapper implements IndividualLogLikelihood {
	
	public WrappedIndividualLogLikelihood(Likelihood originalFunction) {
		super(originalFunction);
	}
	
	@Override
	public Likelihood getOriginalFunction() {return (Likelihood) super.getOriginalFunction();}
	
	@Override
	public Matrix getPredictionVector() {return getOriginalFunction().getPredictionVector();}

	@Override
	public void setYVector(Matrix y) {getOriginalFunction().setYVector(y);}

	@Override
	public Matrix getYVector() {return getOriginalFunction().getYVector();}
	
}
