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

import java.io.Serializable;
import java.security.InvalidParameterException;

import repicea.math.AbstractMathematicalFunctionWrapper;
import repicea.math.MathematicalFunction;
import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.stats.LinearStatisticalExpression;

/**
 * The LinkFunction class is essential to any generalized linear model. It contains a single parameter, which
 * is a LinearStatisticalExpression instance. The two available link functions are the logit and the log link
 * functions.
 * @author Mathieu Fortin - May 2012
 */
@SuppressWarnings("serial")
public class LinkFunction extends AbstractMathematicalFunctionWrapper implements Serializable {
	
	public static enum Type {Logit, Log, CLogLog}
	
	private Type type;
	
	/**
	 * Public constructor. The default nested function is linear.
	 * @param type a Type enum variable that defines the type of link function
	 * @throws InvalidParameterException if type is null
	 */
	public LinkFunction(Type type) {
		this(type, new LinearStatisticalExpression());
	}

	/**
	 * Public constructor.
	 * @param type a Type enum variable that defines the type of link function
	 * @param eta an AbstractMathematicalFunction instance
	 * @throws InvalidParameterException if type is null
	 */
	public LinkFunction(Type type, MathematicalFunction eta) {
		super(eta);
		if (type == null) {
			this.type = Type.Logit;
		}
		this.type = type;
	}

	/**
	 * Provides the type of link function.
	 * @return a Type enum variable
	 */
	public Type getType() {return type;}
	
	@Override
	public Double getValue() {
		double output;
		double parameterEta = getOriginalFunction().getValue();
		switch (type) {
		case Logit:
			output = Math.exp(parameterEta) / (1 + Math.exp(parameterEta));
			return output;
		case Log:
			output = Math.exp(parameterEta);
			return output;
		case CLogLog:
			output = 1 - Math.exp(-Math.exp(parameterEta));
			return output;
		default:
			throw new RuntimeException("This link function has not been implemented : " + type.toString());
		}
	}

	@Override
	public Matrix getGradient() {
		double value;
		double der;
		switch (type) {
		case Logit:
			value = getValue();
			der = value - value * value;
			return getOriginalFunction().getGradient().scalarMultiply(der);
		case Log:
			value = getValue();
			return getOriginalFunction().getGradient().scalarMultiply(value);
		case CLogLog:
			double exp_eta = Math.exp(getOriginalFunction().getValue());
			return getOriginalFunction().getGradient().scalarMultiply(Math.exp(-exp_eta) * exp_eta); 
		default:
			return null;
		}
	}


	@Override
	public SymmetricMatrix getHessian() {
		MathematicalFunction eta = getOriginalFunction();
		double expEta = Math.exp(eta.getValue());
		Matrix gradientProduct = eta.getGradient().multiply(eta.getGradient().transpose());
		
		switch (type) {
		case Logit:
			double denominator = 1d / (1 + expEta);
			double value = getValue();
			double firstDerivative = value - value * value;
			Matrix firstTerm = eta.getHessian().scalarMultiply(firstDerivative);
			
			double secondDerivative = (- expEta * expEta + expEta) * denominator * denominator * denominator;
			Matrix secondTerm = gradientProduct.scalarMultiply(secondDerivative);
		
			return SymmetricMatrix.convertToSymmetricIfPossible(firstTerm.add(secondTerm));
		case Log:
			return SymmetricMatrix.convertToSymmetricIfPossible(gradientProduct.scalarMultiply(expEta).add(eta.getHessian().scalarMultiply(expEta)));
		case CLogLog:
			return SymmetricMatrix.convertToSymmetricIfPossible(gradientProduct.scalarMultiply(-expEta).add(gradientProduct).add(eta.getHessian()).scalarMultiply(Math.exp(-expEta) * expEta));
		default:
			return null;
		}
	}

	
}
