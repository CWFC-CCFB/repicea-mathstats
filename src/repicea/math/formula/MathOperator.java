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

import java.util.HashMap;
import java.util.Map;

/**
 * The MathOperator class is an abstract class for all operators. Basically, an operator has
 * a left side and a right side component as well as a priority. 
 * @author Mathieu Fortin - May 2013
 */
public abstract class MathOperator implements Calculable {

	public static Map<String, Class<? extends MathOperator>> NamedOperators = new HashMap<String, Class<? extends MathOperator>>();
	static {
		NamedOperators.put("exp", Exponential.class);
		NamedOperators.put("log", Logarithm.class);
		NamedOperators.put("sqr", Square.class);
		NamedOperators.put("sqrt", SquareRoot.class);
	}
	
	
	protected int priority;
	protected Calculable leftSide;
	protected Calculable rightSide;

	protected void setLeftSide(Calculable leftSide) {
		this.leftSide = leftSide;
	}
	
	protected void setRightSide(Calculable rightSide) {
		this.rightSide = rightSide;
	}
	
	protected int getPriority() {
		return priority;
	}
	
	/**
	 * This method returns the appropriate mathematical operator.
	 * @param operator a symbol e.g. / * + -
	 * @return a MathOperator instance
	 */
	protected static MathOperator getAppropriateMathOperator(String operator) {
		if (operator.equals("*")) {
			return new MathOperator.Multiply();
		} else if (operator.equals("/")) {
			return new MathOperator.Divide();
		} else if (operator.equals("+")) {
			return new MathOperator.Plus();
		} else if (operator.equals("-")) {
			return new MathOperator.Minus();
		} else if (operator.equals("^")) {
			return new MathOperator.Power();
		} else if (operator.equals("exp")) {
			return new MathOperator.Exponential();
		} else if (operator.equals("log")) {
			return new MathOperator.Logarithm();
		} else if (operator.equals("sqr")) {
			return new MathOperator.Square();
		} else {
			return null;
		}
	}

	protected static String getOperatorLongNameIfAny(Calculable operator) {
		for (String name : NamedOperators.keySet()) {
			if (operator.getClass().isAssignableFrom(NamedOperators.get(name))) {
				return name;
			}
		}
		return "";
	}
	
	/**
	 * A classical plus "+" operator.
	 * @author Mathieu Fortin - May 2013
	 */
	static class Plus extends MathOperator {
		
		protected Plus() {
			priority = 0;
		}
		
 		@Override
		public double calculate() {
			return leftSide.calculate() + rightSide.calculate();
		}
 		
	}

	/**
	 * A classical minus "-" operator.
	 * @author Mathieu Fortin - May 2013
	 */
	static class Minus extends MathOperator {
		
		protected Minus() {
			priority = 0;
		}
		
 		@Override
		public double calculate() {
			return leftSide.calculate() - rightSide.calculate();
		}
 		
	}

	/**
	 * A classical multiply "*" operator.
	 * @author Mathieu Fortin - May 2013
	 */
	static class Multiply extends MathOperator {
		
		protected Multiply() {
			priority = 1;
		}

		@Override
		public double calculate() {
			return leftSide.calculate() * rightSide.calculate();
		}
		
	}

	/**
	 * A classical divide "/" operator.
	 * @author Mathieu Fortin - May 2013
	 */
	static class Divide extends MathOperator {
		
		protected Divide() {
			priority = 1;
		}

		@Override
		public double calculate() {
			return leftSide.calculate() / rightSide.calculate();
		}
		
	}
	
	/**
	 * A classical power "^" operator.
	 * @author Mathieu Fortin - May 2013
	 */
	static class Power extends MathOperator {
		
		protected Power() {
			priority = 2;
		}

		@Override
		public double calculate() {
			return Math.pow(leftSide.calculate(), rightSide.calculate());
		}
	}

	/**
	 * A classical exp operator.
	 * @author Mathieu Fortin - September 2022
	 */
	static class Exponential extends MathOperator {
		
		protected Exponential() {
			priority = 2;
		}

		@Override
		public double calculate() {
			return Math.exp(leftSide.calculate());
		}
	}

	/**
	 * A classical log operator.
	 * @author Mathieu Fortin - September 2022
	 */
	static class Logarithm extends MathOperator {
		
		protected Logarithm() {
			priority = 2;
		}

		@Override
		public double calculate() {
			double argument = leftSide.calculate();
			if (argument <= 0d) {
				throw new UnsupportedOperationException("The MathOperator$Logarithm class has encontered a zero or a negative argument!");
			}
			return Math.log(leftSide.calculate());
		}
	}
	
	/**
	 * A classical square (power 2) operator.
	 * @author Mathieu Fortin - September 2022
	 */
	static class Square extends MathOperator {
		
		protected Square() {
			priority = 2;
		}

		@Override
		public double calculate() {
			double argument = leftSide.calculate();
			return argument * argument;
		}
	}

	/**
	 * A classical square root operator.
	 * @author Mathieu Fortin - May 2023
	 */
	static class SquareRoot extends MathOperator {
		
		protected SquareRoot() {
			priority = 2;
		}

		@Override
		public double calculate() {
			double argument = leftSide.calculate();
			return Math.sqrt(argument);
		}
	}


}
