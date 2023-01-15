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
package repicea.math.formula;

import java.security.InvalidParameterException;
import java.util.Map;

/**
 * The ExpressionHandler class handles either doubles, MathFormula instances or
 * variables.
 * @author Mathieu Fortin - May 2013
 *
 * @param <P> the class that is handled by this instance
 */
abstract class ExpressionHandler<P> implements Calculable {

	protected final P expression;
		
	protected ExpressionHandler(P expression) {
		this.expression = expression;
	}
	
	
	static class DoubleHandler extends ExpressionHandler<Double> {
		
		public DoubleHandler(Double d) {
			super(d);
		}

		@Override
		public double calculate() {
			return expression;
		}
	}
	
	static class MathFormulaHandler extends ExpressionHandler<MathFormula> {

		public MathFormulaHandler(MathFormula expression) {
			super(expression);
		}

		@Override
		public double calculate() {
			return expression.calculate();
		}
		
	}
	
	
	static class VariableHandler extends ExpressionHandler<Map<String, Double>> {

		private final String variable;
		
		protected VariableHandler(Map<String, Double> expression, String variable) {
			super(expression);
			if (!expression.containsKey(variable)) {
				throw new InvalidParameterException("The variable " + variable + " is not part of the function!");
			}
			this.variable = variable;
		}

		@Override
		public double calculate() {
			return expression.get(variable);
		}
		
	}

	static class ParameterHandler extends ExpressionHandler<Map<String, Double>> {

		private final String parameter;
		
		protected ParameterHandler(Map<String, Double> parameterMap, String parameter) {
			super(parameterMap);
			if (!expression.containsKey(parameter)) {
				throw new InvalidParameterException("The parameter " + parameter + " is not part of the function!");
			}
			this.parameter = parameter;
		}

		@Override
		public double calculate() {
			return expression.get(parameter);
		}
		
	}

}
