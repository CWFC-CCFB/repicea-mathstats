/*
 * This file is part of the repicea-statistics library.
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
package repicea.math;

import java.security.InvalidParameterException;
import java.util.Map;

/**
 * The AbstractMathematicalFunctionWrapper class makes it possible to create a function
 * that is based on a nested AbstractMathematicalFunction instance.
 *
 * @author Mathieu Fortin - December 2015
 */
@SuppressWarnings("serial")
public abstract class AbstractMathematicalFunctionWrapper extends AbstractMathematicalFunction {

	private final MathematicalFunction originalFunction;

	/**
	 * Constructor.
	 * @param originalFunction the nested AbstractMathematicalFunction instance
	 */
	public AbstractMathematicalFunctionWrapper(MathematicalFunction originalFunction) {
		if (originalFunction == null) {
			throw new InvalidParameterException("The originalFunction parameter cannot be null");
		}
		this.originalFunction = originalFunction;
	}

	/**
	 * This method returns the original function.
	 * @return an AbstractMathematicalFunction instance
	 */
	public MathematicalFunction getOriginalFunction() {return originalFunction;}
	
	@Override
	public abstract Double getValue();

	@Override
	public abstract Matrix getGradient();

	@Override
	public abstract SymmetricMatrix getHessian();
	

	@Override
	public void setParameterValue(int parameterIndex, double parameterValue) {
		getOriginalFunction().setParameterValue(parameterIndex, parameterValue);
	}

	@Override
	public double getParameterValue(int parameterIndex) {
		return getOriginalFunction().getParameterValue(parameterIndex);
	}

	@Override
	public void setVariableValue(int variableIndex, double variableValue) {
		getOriginalFunction().setVariableValue(variableIndex, variableValue);
	}

	@Override
	public double getVariableValue(int variableIndex) {
		return getOriginalFunction().getVariableValue(variableIndex);
	}

	@Override
	public int getNumberOfParameters() {return getOriginalFunction().getNumberOfParameters();}

	@Override
	public int getNumberOfVariables() {return getOriginalFunction().getNumberOfVariables();}

	
	@Override
	public void setVariables(Matrix xVector) {getOriginalFunction().setVariables(xVector);}
	
	@Override
	public void setParameters(Matrix beta) {getOriginalFunction().setParameters(beta);}
	
	@Override
	public Matrix getParameters() {return getOriginalFunction().getParameters();}

	@Override
	public void setBounds(int parameterIndex, ParameterBound bound) {
		getOriginalFunction().setBounds(parameterIndex, bound);
	}

	@Override
	public boolean isThisParameterValueWithinBounds(int parameterIndex, double parameterValue) {
		return getOriginalFunction().isThisParameterValueWithinBounds(parameterIndex, parameterValue);
	}
}
