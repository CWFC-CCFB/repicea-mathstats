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
package repicea.math;

import java.io.Serializable;
import java.security.InvalidParameterException;
import java.util.HashMap;
import java.util.Map;


/**
 * The AbstractMathematicalFunction class suits mathematical function that can be derived. The return value
 * can be an complex object to enable the nesting of different mathematical functions.
 * @author Mathieu Fortin - June 2011
 */
@SuppressWarnings("serial")
public abstract class AbstractMathematicalFunction implements MathematicalFunction, Serializable {
	
	public static final double MINIMUM_ACCEPTABLE_POSITIVE_VALUE = 1E-8;

	protected final Map<Integer, ParameterBound> parameterBounds;

	private final FastArrayList<Double> parameterValues;
	private final FastArrayList<Double> variableValues;
	
	protected AbstractMathematicalFunction() {
		parameterValues = new FastArrayList<Double>();
		variableValues = new FastArrayList<Double>();
		parameterBounds = new HashMap<Integer, ParameterBound>();
	}

	@Override
	public void setParameterValue(int parameterIndex, double parameterValue) {
		if (parameterBounds.containsKey(parameterIndex)) {
			ParameterBound bound = parameterBounds.get(parameterIndex);
			parameterValue = bound.validateParameter(parameterValue);
		}
		if (parameterIndex < getNumberOfParameters()) {
			parameterValues.set(parameterIndex, parameterValue);
		} else if (parameterIndex == getNumberOfParameters()) {
			parameterValues.add(parameterValue);
		} else {
			throw new InvalidParameterException("The parameter index is not sequential!");
		}
	}

	@Override
	public double getParameterValue(int parameterIndex) {
		return parameterValues.get(parameterIndex);
	}

	@Override
	public void setVariableValue(int variableIndex, double variableValue) {
		if (variableIndex < getNumberOfVariables()) {
			variableValues.set(variableIndex, variableValue);
		} else if (variableIndex == getNumberOfVariables()) {
			variableValues.add(variableValue);
		} else {
			throw new InvalidParameterException("The variable index is not sequential!");
		}
	}
	
	@Override
	public double getVariableValue(int variableIndex) {
		return variableValues.get(variableIndex);
	}
	
	@Override
	public int getNumberOfParameters() {return parameterValues.size();}
	
	@Override
	public int getNumberOfVariables() {return variableValues.size();}

	@Override
	public abstract Double getValue();

	@Override
	public abstract Matrix getGradient();
	
	@Override
	public abstract SymmetricMatrix getHessian();

	@Override
	public void setBounds(int parameterIndex, ParameterBound bound) {
		parameterBounds.put(parameterIndex, bound);
	}

	@Override
	public void setVariables(Matrix xVector) {
		if (!xVector.isRowVector()) {
			throw new IllegalArgumentException("The vector is not a row vector!");
		} else {
			variableValues.clear();
			for (int j = 0; j < xVector.m_iCols; j++) {
				setVariableValue(j, xVector.getValueAt(0, j));
			}
		}
	}
	
	@Override
	public void setParameters(Matrix beta) {
		if (!beta.isColumnVector()) {
			throw new IllegalArgumentException("The vector is not a column vector!");
		} else {
			this.parameterValues.clear();
			for (int i = 0; i < beta.m_iRows; i++) {
				setParameterValue(i, beta.getValueAt(i, 0));
			}
		}
	}
	
	@Override
	public Matrix getParameters() {
		Matrix m = new Matrix(getNumberOfParameters(), 1);
		for (int i = 0; i < getNumberOfParameters(); i++) {
			m.setValueAt(i, 0, getParameterValue(i));
		}
		return m;
	}

	@Override
	public boolean isThisParameterValueWithinBounds(int parameterIndex, double parameterValue) {
		if (parameterBounds.containsKey(parameterIndex)) {
			ParameterBound bound = parameterBounds.get(parameterIndex);
			if (!bound.isParameterValueValid(parameterValue)) {
				return false;
			}
		}
		return true;
	}

//	@Override
//	public Map<Integer, ParameterBound> getBounds() {
//		Map<Integer,ParameterBound> newMap = new HashMap<Integer, ParameterBound>();
//		for (Integer k : parameterBounds.keySet()) {
//			newMap.put(k, parameterBounds.get(k).clone());
//		}
//		return newMap;
//	}
}
