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
package repicea.stats.model.glm.measerr;

import java.security.InvalidParameterException;

import repicea.math.EvaluableFunction;
import repicea.math.MathematicalFunction;
import repicea.math.Matrix;

public abstract class AbstractGLMMeasErrorDefinition implements GLMMeasErrorDefinition {

	static class GradientHessianProvider implements EvaluableFunction<Matrix> {

		private final MathematicalFunction originalFunction;
		private final boolean isForGradient;
		
		GradientHessianProvider(MathematicalFunction originalFunction, boolean isForGradient) {
			this.originalFunction = originalFunction;
			this.isForGradient = isForGradient;
		}
		
		@Override
		public Matrix getValue() {
			return isForGradient ? originalFunction.getGradient() : originalFunction.getHessian();
		}

		@Override
		public void setVariableValue(int variableIndex, double variableValue) {
			originalFunction.setVariableValue(variableIndex, variableValue);
		}

		@Override
		public void setParameterValue(int parameterIndex, double parameterValue) {
			originalFunction.setParameterValue(parameterIndex, parameterValue);
		}

		@Override
		public double getVariableValue(int variableIndex) {
			return originalFunction.getVariableValue(variableIndex);
		}

		@Override
		public double getParameterValue(int parameterIndex) {
			return originalFunction.getParameterValue(parameterIndex);
		}

		@Override
		public void setParameters(Matrix beta) {}

		@Override
		public void setVariables(Matrix xVector) {}

		@Override
		public int getNumberOfParameters() {return originalFunction.getNumberOfParameters();}

		@Override
		public int getNumberOfVariables() {return originalFunction.getNumberOfVariables();}
	}

	protected final MeasurementErrorModel mesErrMod;
	protected final String effectWithMeasError;
	protected int indexEffectWithMeasError;

	/**
	 * Constructor.
	 * @param mesErrMod the type of error model, either Classical or Berkson
	 * @param effectWithMeasError the field name of the variable with measurement error
	 */
	protected AbstractGLMMeasErrorDefinition(MeasurementErrorModel mesErrMod, String effectWithMeasError) {
		this.mesErrMod = mesErrMod;
		if (effectWithMeasError == null || effectWithMeasError.isEmpty()) {
			throw new InvalidParameterException("The effectWithMeasError argument cannot be null or empty!");
		}
		this.effectWithMeasError = effectWithMeasError;
	}

	@Override
	public MeasurementErrorModel getMeasurementErrorModel() {return mesErrMod;}

	public void validate(GLMWithMeasurementError glm) {
		indexEffectWithMeasError = glm.getDataStructure().indexOfThisEffect(effectWithMeasError);
		if (indexEffectWithMeasError == -1) {
			throw new InvalidParameterException("The effect with measurement error " + effectWithMeasError + " is not part of the model definition!");
		} 
	}
}
