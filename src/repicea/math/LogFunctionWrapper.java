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

/**
 * The LogFunctionWrapper class is a specific AbstractMathematicalFunctionWrapper for log functions.
 * 
 * @author Mathieu Fortin - December 2015
 */
@SuppressWarnings("serial")
public class LogFunctionWrapper extends AbstractMathematicalFunctionWrapper {

	/**
	 * Constructor.
	 * @param originalFunction the nested AbstractMathematicalFunction instance
	 */
	public LogFunctionWrapper(MathematicalFunction originalFunction) {
		super(originalFunction);
	}

	@Override
	public Double getValue() {
		return Math.log(getOriginalFunction().getValue());
	}

	@Override
	public Matrix getGradient() {
		return getOriginalFunction().getGradient().scalarMultiply(1d / getOriginalFunction().getValue());
	}

	@Override
	public SymmetricMatrix getHessian() {
		double invValue = 1d/getOriginalFunction().getValue();
		Matrix originalGradient = getOriginalFunction().getGradient();
		Matrix part1 = originalGradient.multiply(originalGradient.transpose()).scalarMultiply(- invValue * invValue);
		SymmetricMatrix part2 = getOriginalFunction().getHessian().scalarMultiply(invValue);
		return SymmetricMatrix.convertToSymmetricIfPossible(part1.add(part2));
	}

}
