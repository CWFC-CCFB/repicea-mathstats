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
 * The LogFunctionWrapper class is a specific AbstractMathematicalFunctionWrapper for exp functions.
 * 
 * @author Mathieu Fortin - December 2015
 */
@SuppressWarnings("serial")
public class ExponentialFunctionWrapper extends	AbstractMathematicalFunctionWrapper {

	public ExponentialFunctionWrapper(MathematicalFunction originalFunction) {
		super(originalFunction);
	}

	@Override
	public Double getValue() {
		return Math.exp(getOriginalFunction().getValue());
	}

	@Override
	public Matrix getGradient() {
		return getOriginalFunction().getGradient().scalarMultiply(getValue());
	}

	@Override
	public SymmetricMatrix getHessian() {
		Matrix part1 = getGradient().multiply(getOriginalFunction().getGradient().transpose());
		SymmetricMatrix part2 = getOriginalFunction().getHessian().scalarMultiply(getValue());
		return SymmetricMatrix.convertToSymmetricIfPossible(part1.add(part2));
	}

}
