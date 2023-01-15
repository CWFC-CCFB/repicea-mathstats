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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;

import repicea.math.formula.ExpressionHandler.DoubleHandler;
import repicea.math.formula.ExpressionHandler.MathFormulaHandler;
import repicea.math.formula.ExpressionHandler.ParameterHandler;
import repicea.math.formula.ExpressionHandler.VariableHandler;

/**
 * The MathFormula class translates a String expression in a tractable mathematical formula.
 * @author Mathieu Fortin - February 2013
 */
public final class MathFormula implements Calculable {
	
	private final static char DOT = new String(".").charAt(0); 
		
	private String formula;
	private final Map<String, Double> parameters;
	private final Map<String, Double> variables;
	
	private final Map<String, MathFormula> nestedMathFormulas;
	private final Map<String, Double> constants;
	private Calculable finalCalculable;
	
	/**
	 * Constructor.
	 * @param formula a String that represents the formula
	 * @param parameters a List of String instances which represent the parameters
	 * @param variables a List of String instances which represent the variables
	 */
	public MathFormula(String formula, LinkedHashMap<String, Double> parameters, LinkedHashMap<String, Double> variables) {
		if (formula == null || formula.isEmpty()) {
			throw new InvalidParameterException("The formula cannot be null or empty!");
		}
		
		nestedMathFormulas = new HashMap<String, MathFormula>();
		constants = new HashMap<String, Double>();
		
		this.formula = formula.replace(" ", "").trim();
		this.parameters = new LinkedHashMap<String, Double>();
		if (parameters != null) {
			this.parameters.putAll(parameters);
		}
		this.variables = new LinkedHashMap<String, Double>();
		if (variables != null) {
			this.variables.putAll(variables);
		}
		parseParentheses();
		checkParametersAndVariables();
		defineOperators();
	}

	/**
	 * Constructor for nested formulas.
	 * @param formula a String that represents the formula
	 * @param operator a long named operator if needed (can be null if none)
	 * @param parentFormula the MathFormula instance that is creating this nested formula
	 */
	private MathFormula(String formula, MathOperator operator, MathFormula parentFormula) {
		if (formula == null || formula.isEmpty()) {
			throw new InvalidParameterException("The formula argument should be non null and non empty!");
		}
		
		nestedMathFormulas = new HashMap<String, MathFormula>();
		constants = new HashMap<String, Double>();
		
		this.formula = formula.replace(" ", "").trim();
		this.parameters = parentFormula.parameters;
		this.variables = parentFormula.variables;
		
		parseParentheses();
		checkParametersAndVariables();
		this.finalCalculable = operator;
		defineOperators();
	}


	private void defineOperators() {
		String formula = this.formula;
	
		StringTokenizer tkz = new StringTokenizer(formula, "*/+-^");
		
		String leftPart = tkz.nextToken().trim();
		ExpressionHandler<?> leftExpression = getAppropriateExpressionHandler(leftPart);
		
		String rightPart;
		ExpressionHandler<?> rightExpression = null;
		MathOperator currentMathOperator;
		MathOperator lastMathOperator = null;
		
		while (tkz.hasMoreTokens()) {
			
			rightPart = tkz.nextToken().trim();
			rightExpression = getAppropriateExpressionHandler(rightPart);
			
			int index0 = formula.indexOf(leftPart) + leftPart.length();
			int index1 = formula.indexOf(rightPart);
			String currentOperator = formula.substring(index0, index1).trim();
			
			formula = formula.substring(formula.indexOf(currentOperator) + 1);
			currentMathOperator = MathOperator.getAppropriateMathOperator(currentOperator);
			
			if (lastMathOperator != null) {
				if (lastMathOperator.getPriority() >= currentMathOperator.getPriority()) {
					currentMathOperator.setLeftSide(lastMathOperator);
					currentMathOperator.setRightSide(rightExpression);
					finalCalculable = currentMathOperator;
				} else {
					lastMathOperator.setRightSide(currentMathOperator);
					currentMathOperator.setLeftSide(leftExpression);
					currentMathOperator.setRightSide(rightExpression);
				}
			} else {
				if (finalCalculable != null) {
					currentMathOperator.setLeftSide(leftExpression);
					currentMathOperator.setRightSide(rightExpression);
					((MathOperator) finalCalculable).setLeftSide(currentMathOperator);
				} else {
					finalCalculable = currentMathOperator;
					currentMathOperator.setLeftSide(leftExpression);
					currentMathOperator.setRightSide(rightExpression);
				}
			}
			lastMathOperator = currentMathOperator;
			leftPart = rightPart;
			leftExpression = rightExpression;
		}
		
		if (finalCalculable == null) {		// it means there was no operator at all in the formula e.g. ($nf3)
			finalCalculable = leftExpression;
		} else if (rightExpression == null) {
			((MathOperator) finalCalculable).setLeftSide(leftExpression);
		}

	}
	
	private ExpressionHandler<?> getAppropriateExpressionHandler(String expression) {
		if (constants.containsKey(expression)) {
			return new DoubleHandler(constants.get(expression));
		} else if (nestedMathFormulas.containsKey(expression)) {
			return new MathFormulaHandler(nestedMathFormulas.get(expression));
		} else if (parameters.containsKey(expression)) {
			return new ParameterHandler(parameters, expression);		
		} else if (variables.containsKey(expression)) {
			return new VariableHandler(variables, expression);
		} else {
			return null;
		}
	}
	
	private void replaceInFormula(String str, String id) {
		int startIndex = formula.indexOf(str);
		int endIndex = startIndex + str.length();
		String newFormula = formula.substring(0,startIndex) + id + formula.substring(endIndex);
		formula = newFormula;
	}
	
	
	/**
	 * This method makes nested math formulas from expressions in parentheses.
	 */
	private void parseParentheses() {
		int firstOccurrence;
		int firstOccurrenceClosing;
		while ((firstOccurrence = formula.indexOf("(")) != -1) {
			firstOccurrenceClosing = formula.indexOf(")");
			if (firstOccurrenceClosing < firstOccurrence) {
				throw new InvalidParameterException("The closing parenthesis comes before the opening one in the formula!");
			} 
			int numberOfOpeningParentheses = 0;
			Character charact;
			for (int i = firstOccurrence + 1; i < formula.length(); i++) {
				charact = formula.charAt(i);
				if (charact.toString().equals("(")) {
					numberOfOpeningParentheses++;
				}
				if (charact.toString().equals(")")) {
					if (numberOfOpeningParentheses == 0) {
						String nestedFormula = formula.substring(firstOccurrence + 1, i);
						String longNamedOperator = checkIfLongNamedOperatorBefore(formula, firstOccurrence);
						
						MathFormula nestedMathFormula = longNamedOperator != null ?
								new MathFormula(nestedFormula, MathOperator.getAppropriateMathOperator(longNamedOperator), this) :
								new MathFormula(nestedFormula, null, this);

						String nestedFormulaID = "$nf" + nestedMathFormulas.size();
						nestedMathFormulas.put(nestedFormulaID, nestedMathFormula);
						int firstCut = longNamedOperator != null ?
								firstOccurrence - longNamedOperator.length() :
									firstOccurrence;
						String newFormula = formula.substring(0, firstCut) + nestedFormulaID + formula.substring(i + 1);
						formula = newFormula;
						break;
					} else {
						numberOfOpeningParentheses--;
					}
				}
			}
			if (numberOfOpeningParentheses != 0) {
				throw new InvalidParameterException("The number of closing parentheses does not match the number of opening parentheses!");
			}
		}
		firstOccurrenceClosing = formula.indexOf(")");
		if (firstOccurrenceClosing != -1) {	
			throw new InvalidParameterException("Either the opening or closing parenthesis is missing!");
		} 
	}


	
	private String checkIfLongNamedOperatorBefore(String thisFormula, int firstOccurrence) {
		for (String longNamedOperator : MathOperator.NamedOperators.keySet()) {
			if (firstOccurrence - longNamedOperator.length() > -1) {
				String thatPossibleOperator = thisFormula.substring(firstOccurrence - longNamedOperator.length(), firstOccurrence);
				if (longNamedOperator.equals(thatPossibleOperator)) {
					return longNamedOperator;
				}
			}
		}
		return null;
	}


	@Override
	public String toString() {
		String originalFormula = formula;
		MathFormula nestedFormula;
		for (String key : nestedMathFormulas.keySet()) {
			nestedFormula = nestedMathFormulas.get(key);
			String longOperatorIfAny = MathOperator.getOperatorLongNameIfAny(nestedFormula.finalCalculable);
			originalFormula = originalFormula.replace(key, longOperatorIfAny +  "(" + nestedFormula.toString() + ")");
		}
		for (String key : constants.keySet()) {
			double d = constants.get(key);
			originalFormula = originalFormula.replace(key, ((Double) d).toString());
		}
		return originalFormula;
	}
	
	private boolean isTokenKnown(String token) {
		if (token.startsWith("$"))
			return true;
		for (String longNamedOperator : MathOperator.NamedOperators.keySet()) {
			if (token.startsWith(longNamedOperator + "$"))
				return true;
		}
		return false;
	}

	private void checkParametersAndVariables() {
		List<String> constantList = new ArrayList<String>();
		StringTokenizer tkz = new StringTokenizer(formula, "+-*/^");
		String token;
		while (tkz.hasMoreTokens()) {
			token = tkz.nextToken().trim();
			if (!isTokenKnown(token)) {
				if (isNumeric(token)) {
					constantList.add(token);
				} else {
					if (!parameters.containsKey(token) && !variables.containsKey(token)) {
						throw new InvalidParameterException("The expression " + token + " has not been defined as parameter or variable!");
					}
				}
			}
		}
		
		for (String constant : constantList) {
			registerConstant(constant);
		}
		
	}
	
	private void registerConstant(String token) {
		String constantID = "$c" + constants.size();
		constants.put(constantID, Double.parseDouble(token));
		replaceInFormula(token, constantID);
	}

	private boolean isNumeric(String str) {
		boolean hasDot = false;
		char c;
		for (int i = 0; i < str.length(); i++) {
			c = str.charAt(i);
			if (!Character.isDigit(c)) {
				if (DOT != c) {
					return false;
				} else {
					if (hasDot) {
						return false;
					} else {
						hasDot = true;
					}
				
				}
			}
		}
		return true;
	}
	

	@Override
	public double calculate() {
		return finalCalculable.calculate();
	}

	public void setParameter(String name, double value) {
		if (!parameters.containsKey(name) ) {
			throw new InvalidParameterException("The parameter " + name + " has not been defined!");
		}
		parameters.put(name, value);
	}

	public void setVariable(String name, double value) {
		if (!variables.containsKey(name) ) {
			throw new InvalidParameterException("The variable " + name + " has not been defined!");
		}
		variables.put(name, value);
	}

	public List<String> getVariables() {
		return new ArrayList<String>(variables.keySet());
	}


}
