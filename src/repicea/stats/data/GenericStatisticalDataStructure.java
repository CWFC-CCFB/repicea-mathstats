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

	
	@SuppressWarnings({ "rawtypes", "unchecked" })
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
			StringTokenizer tkzExclusiveInteraction = new StringTokenizer(effectName, ":");
			effectsInThisInteraction.clear();
			while (tkzExclusiveInteraction.hasMoreTokens()) {
				effectsInThisInteraction.add(tkzExclusiveInteraction.nextToken());
			}
			Matrix matXForThisEffect = null;
			for (String singleEffect : effectsInThisInteraction) {
				Matrix matXtmp;
				MathFormula formula;
				if ((formula = effects.get(singleEffect)) != null) {	// there is a Math Formula behind this effect
					String fName = formula.getVariables().get(0);
					Matrix originalValues = getVectorOfThisField(fName);
					matXtmp = new Matrix(originalValues.m_iRows, 1);
					for (int i = 0; i < originalValues.m_iRows; i++) {
						formula.setVariable(fName, originalValues.getValueAt(i, 0));
						matXtmp.setValueAt(i, 0, formula.calculate());
					}
				} else {
					int indexOfReferenceClass = singleEffect.indexOf("#");
					String refClass = null;
					if (indexOfReferenceClass != -1) {
						refClass = singleEffect.substring(indexOfReferenceClass + 1);
						singleEffect = singleEffect.substring(0, indexOfReferenceClass);
					}

					Class<?> fieldType = dataSet.getFieldTypeOfThisField(singleEffect);
					if (Number.class.isAssignableFrom(fieldType)) {		// it is either a double or an integer
						matXtmp = getVectorOfThisField(singleEffect);
					} else {
						matXtmp = computeDummyVariables(singleEffect, refClass);
					}
				}
				matXForThisEffect = matXForThisEffect == null ?
						matXtmp :
							MatrixUtility.combineMatrices(matXForThisEffect, matXtmp);
			}
			matrixX = matrixX == null ? 
					matXForThisEffect :
						matrixX.matrixStack(matXForThisEffect, false);
		}
		return matrixX;
	}

	@Override
	public Matrix constructVectorY() {return dataSet.getVectorOfThisField(yName);}
	
	private List<String> getLongNamedEffects(String modelEffects) {
		List<String> longNamedEffects = new ArrayList<String>();
		for (String longNamedOperator : MathOperator.NamedOperators.keySet()) {
			List<String> longNamedEffectsForThisOperator = ObjectUtility.extractSequences(modelEffects, longNamedOperator + "(", ")");
			if (longNamedEffectsForThisOperator.size() > 1) {
				for (int i = 1; i < longNamedEffectsForThisOperator.size(); i++)
					longNamedEffects.add(longNamedEffectsForThisOperator.get(i));
			}
		}
		return longNamedEffects;
	}
	
	@Override
	public void setModelDefinition(String modelDefinition) {
		List<String> responseAndFixedEffects = ObjectUtility.decomposeUsingToken(modelDefinition, "~"); 
		
		if (responseAndFixedEffects.size() != 2) {
			throw new InvalidParameterException("The model specification is incorrect!");
		}

		yName = responseAndFixedEffects.get(0);
		String modelEffects = responseAndFixedEffects.get(1);

//		List<String> longNamedEffects = new ArrayList<String>();
//		for (String longNamedOperator : MathOperator.NamedOperators.keySet()) {
//			List<String> longNamedEffectsForThisOperator = ObjectUtility.extractSequences(modelEffects, longNamedOperator + "(", ")");
//			if (longNamedEffectsForThisOperator.size() > 1) {
//				for (int i = 1; i < longNamedEffectsForThisOperator.size(); i++)
//					longNamedEffects.add(longNamedEffectsForThisOperator.get(i));
//			}
//		}
		List<String> longNamedEffects = getLongNamedEffects(modelEffects);

		Map<String, String> longNamedEffectsMap = new HashMap<String, String>();
		int id = 0;
		String substitute;
		for (String longNamedEffect : longNamedEffects) {
			substitute = "&" + id;
			modelEffects = modelEffects.replace(longNamedEffect, substitute);
			longNamedEffectsMap.put(substitute, longNamedEffect);
			id++;
		}
		
		effects.clear();
		List<String> effectsInThisInteraction = new ArrayList<String>();
		List<String> effectAndInteractionList = ObjectUtility.decomposeUsingToken(modelEffects, "+");

		for (String effectOrInteraction : effectAndInteractionList) {
			StringTokenizer tkzInclusiveInteraction = new StringTokenizer(effectOrInteraction, "*");
			StringTokenizer tkzExclusiveInteraction = new StringTokenizer(effectOrInteraction, ":");
			effectsInThisInteraction.clear();

			int numberOfInclusiveInteraction = tkzInclusiveInteraction.countTokens();
			int numberOfExclusiveInteraction = tkzExclusiveInteraction.countTokens();

			if (numberOfInclusiveInteraction > 1 && numberOfExclusiveInteraction > 1) {
				throw new StatisticalDataException("Error : symbols * and : are being used at the same time in the model specification!");
			}

			boolean isAnInclusiveInteraction = numberOfInclusiveInteraction > 1;
			if (isAnInclusiveInteraction) {
				while (tkzInclusiveInteraction.hasMoreTokens()) {
					effectsInThisInteraction.add(tkzInclusiveInteraction.nextToken());
				}
				Collections.sort(effectsInThisInteraction);
				effectsInThisInteraction.addAll(getAllCombinations(effectsInThisInteraction));
			} else {
				List<String> tempColl = new ArrayList<String>();
				while (tkzExclusiveInteraction.hasMoreTokens()) {
					tempColl.add(tkzExclusiveInteraction.nextToken());
				}
				Collections.sort(tempColl);
				String interaction = "";
				for (int i = 0; i < tempColl.size(); i++)
					interaction = i == tempColl.size() - 1 ? 
							interaction + tempColl.get(i) :
								interaction + tempColl.get(i) + ":";
				effectsInThisInteraction.add(interaction);
			}
			for (String effectName : effectsInThisInteraction) {
				if (!effects.containsKey(effectName)) {
					MathFormula formula = null;
					if (longNamedEffectsMap.containsKey(effectName)) {
						String originalEffectName = longNamedEffectsMap.get(effectName);
						formula = extractFormulaIfAny(originalEffectName);
					}
					effects.put(effectName, formula);
				}
			}
		}
	}
	
	protected static final List<String> getAllCombinations(List<String> simpleEffects) {
		List<String> outputList = new ArrayList<String>();
		List<String> uniqueValues = new ArrayList<String>();
		for (String str : simpleEffects) {
			if (!uniqueValues.contains(str)) {
				uniqueValues.add(str);
			}
		}
		if (uniqueValues.size() == 2) {
			outputList.add(uniqueValues.get(0) + ":" + uniqueValues.get(1)); 
		} else if (uniqueValues.size() == 3) {
			for (int i = 0; i < uniqueValues.size() - 1; i++) {
				for (int j = i+1; j < uniqueValues.size(); j++) {
					outputList.add(uniqueValues.get(i) + ":" + uniqueValues.get(j));
				}
			}
			outputList.add(uniqueValues.get(0) + ":" + uniqueValues.get(1) + ":" + uniqueValues.get(2));
		} else {
			throw new InvalidParameterException("Interactions with more than 3 effects are not allowed!");
		}
		return outputList;
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
