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
package repicea.math;

import java.security.InvalidParameterException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.TreeMap;

/**
 * The InternalMathematicalFunctionWrapper class wraps a MathematicalFunction instance and 
 * redefines the parameter and variable indices. <br>
 * <br>
 * It allows for product of functions for instance.
 * @author Mathieu Fortin - July 2022
 */
@SuppressWarnings("serial")
public class InternalMathematicalFunctionWrapper extends AbstractMathematicalFunctionWrapper {

	final TreeMap<Integer, Integer> parmMap;		// from new to original
	final TreeMap<Integer, Integer> reverseParmMap; // from original to new
	final TreeMap<Integer, Integer> varMap;			// from new to original
	final TreeMap<Integer, Integer> reverseVarMap;  // from original to new
	
	/**
	 * Constructor.
	 * @param originalFunction the MathematicalFunction instance to be wrapped
	 * @param newParmIndex the new parameter indices 
	 * @param newVariableIndex the new variable indices
	 */
	public InternalMathematicalFunctionWrapper(MathematicalFunction originalFunction, List<Integer> newParmIndex, List<Integer> newVariableIndex) {
		super(originalFunction);
		parmMap = new TreeMap<Integer, Integer>();
		varMap = new TreeMap<Integer, Integer>();
		reverseParmMap = new TreeMap<Integer, Integer>();
		reverseVarMap = new TreeMap<Integer, Integer>();
		if (newParmIndex == null || newParmIndex.size() != getOriginalFunction().getNumberOfParameters()) {
			throw new InvalidParameterException("The size of the newParmIndex argument is not consistent with the number of parameters in the original function!");
		}
		if (newVariableIndex == null || newVariableIndex.size() != getOriginalFunction().getNumberOfVariables()) {
			throw new InvalidParameterException("The size of the newVariableIndex argument is not consistent with the number of variables in the original function!");
		}
		for (Integer index : newParmIndex) {
			parmMap.put(index, newParmIndex.indexOf(index));
			reverseParmMap.put(newParmIndex.indexOf(index), index);
		}
		for (Integer index : newVariableIndex) {
			varMap.put(index, newVariableIndex.indexOf(index));
			reverseVarMap.put(newVariableIndex.indexOf(index), index);
		}
	}

	@Override
	public final void setParameterValue(int index, double value) {
		if (!parmMap.containsKey(index)) {
			throw new InvalidParameterException("The parameter index " + index + " is not valid!");
		} else {
			super.setParameterValue(parmMap.get(index), value);
		}
	}

	@Override
	public final double getParameterValue(int index) {
		if (!parmMap.containsKey(index)) {
			throw new InvalidParameterException("The parameter index " + index + " is not valid!");
		} else {
			return super.getParameterValue(parmMap.get(index));
		}
	}
	
	@Override
	public void setVariableValue(int index, double value) {
		if (!varMap.containsKey(index)) {
			throw new InvalidParameterException("The variable index " + index + " is not valid!");
		} else {
			super.setVariableValue(varMap.get(index), value);
		}
	}
	
	@Override
	public final double getVariableValue(int index) {
		if (!varMap.containsKey(index)) {
			throw new InvalidParameterException("The variable index " + index + " is not valid!");
		} else {
			return super.getVariableValue(varMap.get(index));
		}
	}

	@Override
	public Double getValue() {return getOriginalFunction().getValue();}

	@Override
	public Matrix getGradient() {return getOriginalFunction().getGradient();}

	@Override
	public SymmetricMatrix getHessian() {return getOriginalFunction().getHessian();}

	Collection<Integer> getNewParameterIndices() {
		return parmMap.keySet();
	}
	
	Collection<Integer> getNewVariableIndices() {
		return varMap.keySet();
	}

	/**
	 * Produce a list of integers. 
	 * 
	 * @param from the initial value of the list.
	 * @param to the final value of the list.
	 * @return a List instance with Integers
	 */
	public static List<Integer> produceListFromTo(int from, int to) {
		if (to < from) {
			throw new InvalidParameterException("The to argument should be equal to or greater than the from argument!");
		}
		List<Integer> myList = new ArrayList<Integer>();
		for (int i = from; i <= to; i++) {
			myList.add(i);
		}
		return myList;
	}

	@Override
	public boolean isThisParameterValueWithinBounds(int parameterIndex, double parameterValue) {
		if (!parmMap.containsKey(parameterIndex)) {
			throw new InvalidParameterException("This parameter index is invalid!");
		} else {
			return getOriginalFunction().isThisParameterValueWithinBounds(parmMap.get(parameterIndex), parameterValue);
		}
	}
	
	@Override 
	public void setBounds(int parameterIndex, ParameterBound bounds) {
		if (!parmMap.containsKey(parameterIndex)) {
			throw new InvalidParameterException("This parameter index is invalid!");
		} else {
			getOriginalFunction().setBounds(parmMap.get(parameterIndex), bounds);
		}
	}
	
}
