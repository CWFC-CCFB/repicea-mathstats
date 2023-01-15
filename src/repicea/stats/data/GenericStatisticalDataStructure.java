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
package repicea.stats.data;

import java.security.InvalidParameterException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;
import java.util.Vector;

import repicea.math.Matrix;
import repicea.math.formula.MathFormula;
import repicea.math.formula.MathOperator;
import repicea.math.utility.MatrixUtility;
import repicea.util.ObjectUtility;

/**
 * The StatisticalDataStructure class is an abstract class that implements all the features to be able
 * to fit a statistical model. The structure includes a vector of dependent variables (vectorY), a matrix
 * of covariates (matrixX) and a possible hierarchical data structure (hierarchicalStructure).
 * @author Mathieu Fortin - June 2011
 */
public class GenericStatisticalDataStructure implements StatisticalDataStructure {

	protected final DataSet dataSet;
	protected boolean isInterceptModel;
//	protected Matrix vectorY;
//	protected Matrix matrixX;
	protected final LinkedHashMap<String, MathFormula> effects;
	protected String yName;
	
	/**
	 * General constructor.
	 * @param dataSet the DataSet instance from which the structure is going to be extracted
	 */
	public GenericStatisticalDataStructure(DataSet dataSet) {
		this.dataSet= dataSet;
		isInterceptModel = true;
		effects = new LinkedHashMap<String, MathFormula>();
	}
	

	@SuppressWarnings({ "rawtypes"})
	protected Matrix computeDummyVariables(String fieldName, String refClass) throws StatisticalDataException {
		List possibleValues = getPossibleValueForDummyVariable(fieldName, refClass);
		int fieldIndex = getDataSet().getIndexOfThisField(fieldName);
		return dataSet.getDummyMatrix(possibleValues, fieldIndex);
	}

	
	@Override
	public List getPossibleValueForDummyVariable(String fieldName, String refClass) {
		int fieldIndex = getDataSet().getIndexOfThisField(fieldName);
		if (fieldIndex ==  -1) {
			throw new InvalidParameterException("Field " + fieldName + " is not part of the DataSet instance!");
		}
		List possibleValues = dataSet.getPossibleValuesInThisField(fieldIndex);
		Collections.sort(possibleValues);
		if (refClass != null && !possibleValues.contains(refClass)) {
			throw new InvalidParameterException("Reference class category " + refClass + " does not belong to this class variable!");
		}
				
		if (isInterceptModel()) {
			if (refClass != null) {
				possibleValues.remove(possibleValues.indexOf(refClass));
			} else {
				possibleValues.remove(0);
			}
		} 
		return possibleValues;
	}
	
	
	@Override
	public boolean isInterceptModel() {return isInterceptModel;}
	
	@Override
	public void setInterceptModel(boolean isInterceptModel) {this.isInterceptModel = isInterceptModel;}
	
	protected Matrix getVectorOfThisField(String fName) {
		return dataSet.getVectorOfThisField(fName);
	}
	
	@Override
	public Matrix constructMatrixX() {
		Matrix matrixX = null;
		
		if (this.isInterceptModel) {
			matrixX = new Matrix(getNumberOfObservations(), 1, 1, 0);
		}
		
		Vector<String> effectsInThisInteraction = new Vector<String>();

		for (String effectName : effects.keySet()) {
			Matrix subMatrixX = null;
			MathFormula formula;
			if ((formula = effects.get(effectName)) != null) {	// there is a Math Formula behind this effect
				String fName = formula.getVariables().get(0);
				Matrix originalValues = getVectorOfThisField(fName);
				subMatrixX = new Matrix(originalValues.m_iRows, 1);
				for (int i = 0; i < originalValues.m_iRows; i++) {
					formula.setVariable(fName, originalValues.getValueAt(i, 0));
					subMatrixX.setValueAt(i, 0, formula.calculate());
				}
			} else {
				StringTokenizer tkzInclusiveInteraction = new StringTokenizer(effectName, "*");
				StringTokenizer tkzExclusiveInteraction = new StringTokenizer(effectName, ":");
				effectsInThisInteraction.clear();

				int numberOfInclusiveInteraction = tkzInclusiveInteraction.countTokens();
				int numberOfExclusiveInteraction = tkzExclusiveInteraction.countTokens();

				if (numberOfInclusiveInteraction > 1 && numberOfExclusiveInteraction > 1) {
					throw new StatisticalDataException("Error : symbols * and : are being used at the same time in the model specification!");
				}

				boolean isAnInclusiveInteraction = numberOfInclusiveInteraction > 1;
				StringTokenizer selectedTokenizer;
				if (isAnInclusiveInteraction) {
					selectedTokenizer = tkzInclusiveInteraction;
				} else {
					selectedTokenizer = tkzExclusiveInteraction;
				}

				while (selectedTokenizer.hasMoreTokens()) {
					effectsInThisInteraction.add(selectedTokenizer.nextToken());
				}

				Matrix matrixTmp;
				for (String effect : effectsInThisInteraction) {

					int indexOfReferenceClass = effect.indexOf("#");
					String refClass = null;
					if (indexOfReferenceClass != -1) {
						refClass = effect.substring(indexOfReferenceClass + 1);
						effect = effect.substring(0, indexOfReferenceClass);
					}

					Class<?> fieldType = dataSet.getFieldTypeOfThisField(effect);
					if (Number.class.isAssignableFrom(fieldType)) {		// it is either a double or an integer
						matrixTmp = getVectorOfThisField(effect);
					} else {
						matrixTmp = computeDummyVariables(effect, refClass);
					}

					if (subMatrixX == null) {
						subMatrixX = matrixTmp;
					} else {
						subMatrixX = MatrixUtility.combineMatrices(subMatrixX, matrixTmp);
					}
				}
			}

			if (matrixX == null) {
				matrixX = subMatrixX;
			} else {
				matrixX = matrixX.matrixStack(subMatrixX, false);
			}
		}
		return matrixX;
	}

	@Override
	public Matrix constructVectorY() {return dataSet.getVectorOfThisField(yName);}
	
	@Override
	public void setModelDefinition(String modelDefinition) {
		List<String> responseAndFixedEffects = ObjectUtility.decomposeUsingToken(modelDefinition, "~"); 
		
		if (responseAndFixedEffects.size() != 2) {
			throw new InvalidParameterException("The model specification is incorrect!");
		}

		yName = responseAndFixedEffects.get(0);
		String modelEffects = responseAndFixedEffects.get(1);

		List<String> longNamedEffects = new ArrayList<String>();
		for (String longNamedOperator : MathOperator.NamedOperators.keySet()) {
			List<String> longNamedEffectsForThisOperator = ObjectUtility.extractSequences(modelEffects, longNamedOperator + "(", ")");
			if (longNamedEffectsForThisOperator.size() > 1) {
				for (int i = 1; i < longNamedEffectsForThisOperator.size(); i++)
					longNamedEffects.add(longNamedEffectsForThisOperator.get(i));
			}
		}

		Map<String, String> longNamedEffectsMap = new HashMap<String, String>();
		int id = 0;
		String substitute;
		for (String longNamedEffect : longNamedEffects) {
			substitute = "&" + id;
			modelEffects = modelEffects.replace(longNamedEffect, substitute);
			longNamedEffectsMap.put(substitute, longNamedEffect);
			id++;
		}
		
		List<String> effectList = ObjectUtility.decomposeUsingToken(modelEffects, "+");
		effects.clear();

		for (String effectName : effectList) {
			MathFormula formula = null;
			if (longNamedEffectsMap.containsKey(effectName)) {
				String originalEffectName = longNamedEffectsMap.get(effectName);
				formula = extractFormulaIfAny(originalEffectName);
			}
			effects.put(effectName, formula);
		}
	}
	
	private MathFormula extractFormulaIfAny(String effectName) {
		boolean isLongNamedOperator = false;
		for (String longNamedOperator : MathOperator.NamedOperators.keySet()) {
			if (effectName.startsWith(longNamedOperator + "(")) {
				isLongNamedOperator = true;
				break;
			}
		}
		if (isLongNamedOperator) {
			LinkedHashMap<String, Double> variables = new LinkedHashMap<String, Double>();
			for (String fieldName : dataSet.getFieldNames()) {
				if (effectName.contains(fieldName)) {
					variables.put(fieldName, 0d);
				}
			}
			return new MathFormula(effectName, null, variables);
		} else {
			return null;
		}
	}


	@Override
	public int getNumberOfObservations() {return dataSet.getNumberOfObservations();}


	@Override
	public DataSet getDataSet() {return dataSet;}

	@Override
	public int indexOfThisEffect(String effect) {
		int index = getEffectList().indexOf(effect);
		if (index != -1 && isInterceptModel()) {
			return index + 1;
		} else {
			return index;
		}
	}
	
	@Override
	public List<String> getEffectList() {
		List<String> effectList = new ArrayList<String>();
		for (String key : effects.keySet()) {
			MathFormula f = effects.get(key);
			effectList.add(f != null ? f.toString() : key);
		}
		return effectList;
	}
}
