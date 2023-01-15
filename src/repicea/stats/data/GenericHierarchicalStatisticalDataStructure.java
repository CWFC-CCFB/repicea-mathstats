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
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import repicea.math.Matrix;
import repicea.util.ObjectUtility;

/**
 * This class is the basic class for hierarchical data structure. 
 * @author Mathieu Fortin - October 2011
 */
public class GenericHierarchicalStatisticalDataStructure extends GenericStatisticalDataStructure implements HierarchicalStatisticalDataStructure {
	
	protected final Map<String, DataBlock> hierarchicalStructure;		// from outer to inner levels
	protected Map<String, List<String>> randomEffectsSpecifications;

	
	/**
	 * The string represents the hierarchical level.
	 */
	protected final Map<String, Matrix> matricesZ;

	protected final boolean sorted;
	
	/**
	 * General constructor. To be defined in derived class.
	 * @param dataSet a DataSet instance
	 * @param sorted a boolean. If true the dataset will be sorted according to the hierarchical structure
	 */
	public GenericHierarchicalStatisticalDataStructure(DataSet dataSet, boolean sorted) {
		super(dataSet);
		this.sorted = sorted;
		hierarchicalStructure = new LinkedHashMap<String, DataBlock>();
		matricesZ = new LinkedHashMap<String, Matrix>();
	}

	/**
	 * Constructor for hierarchically sorted dataset.
	 * @param dataSet a DataSet instance
	 */
	public GenericHierarchicalStatisticalDataStructure(DataSet dataSet) {
		this(dataSet, true);
	}
	
	@Override
	public void setModelDefinition(String modelDefinition) {
		modelDefinition = extractRandomEffects(modelDefinition);
		super.setModelDefinition(modelDefinition);
	}


	/**
	 * This method filters the model definition for sub strings in parentheses. These substrings are then
	 * recorded as random effects and taken out of the original model definition.
	 * @param modelDefinition the model definition
	 * @return the model definition without the random effects
	 * @throws StatisticalDataException 
	 */
	private String extractRandomEffects(String modelDefinition) {
		List<String> occurrences = ObjectUtility.extractSequences(modelDefinition, "(", ")");
		List<String> retainedOccurrences = new ArrayList<String>();
		for (String possOcc : occurrences) {
			if (possOcc.contains("|")) {
				retainedOccurrences.add(possOcc);
			}
		}
		if (retainedOccurrences.size() > 1) {
			throw new InvalidParameterException("The model definition only supports one random effect statement!");
		}
		
		if (!retainedOccurrences.isEmpty()) {
			for (String randomEffect : retainedOccurrences) {
				recordRandomEffects(randomEffect);
			}
			
			List<String> hierarchicalLevels = new ArrayList<String>();
			for (String level : getRandomEffectsSpecifications().keySet()) {
				hierarchicalLevels.add(level);
			}
			
			if (sorted) {
				sortDataAccordingToRandomEffects(hierarchicalLevels);		// we sort the data before setting the hierarchical structure
			}
			
			setHierarchicalStructureLevel(hierarchicalLevels);

			setRandomEffectStructure();
		}

		
		return retainedOccurrences.isEmpty() ? // means there is no random effect specification 
				modelDefinition : 
					modelDefinition.replace(retainedOccurrences.get(0), "");
	}
	
		
	private void sortDataAccordingToRandomEffects(List<String> hierarchicalLevels) {
		List<Integer> indexList = new ArrayList<Integer>();
		for (String levelName : hierarchicalLevels) {
			indexList.add(dataSet.getIndexOfThisField(levelName));
		}
		dataSet.sortObservations(indexList);
	}


	protected void recordRandomEffects(String effectName) {
		String randomEffectSpec = effectName.replace("(","").replace(")","");
		List<String> randomEffectComponents = ObjectUtility.decomposeUsingToken(randomEffectSpec, "|");
		if (randomEffectComponents.size() != 2) {
			throw new InvalidParameterException("The random effect " + effectName + " is not properly defined!");
		}
		
		List<String> hierarchicalLevels = ObjectUtility.decomposeUsingToken(randomEffectComponents.get(1), "/");
		for (String level : hierarchicalLevels) {
			if (dataSet.getIndexOfThisField(level) == -1) {
				throw new InvalidParameterException("Field " + level + "does not exist!");
			}
		}
		
		List<String> randomEffectsForTheseLevels = ObjectUtility.decomposeUsingToken(randomEffectComponents.get(0), "+");
		for (String randomEffect : randomEffectsForTheseLevels) {
			if (!randomEffect.equals("1")) {
				if (dataSet.getIndexOfThisField(randomEffect) == -1) {
					throw new InvalidParameterException("Field " + randomEffect + "does not exist!");
				}
			}
		}
		
		List<String> effects;
		for (String level : hierarchicalLevels) {
			if (!getRandomEffectsSpecifications().containsKey(level)) {
				getRandomEffectsSpecifications().put(level, new ArrayList<String>());
			} 
			effects = getRandomEffectsSpecifications().get(level);
			for (String effect : randomEffectsForTheseLevels) {
				if (!effects.contains(effect)) {
					effects.add(effect);
				}
			}
		}
	}
	
	private Map<String, List<String>> getRandomEffectsSpecifications() {
		if (randomEffectsSpecifications == null) {
			randomEffectsSpecifications = new LinkedHashMap<String, List<String>>();
		}
		return randomEffectsSpecifications;
	}
	
	protected void setRandomEffectStructure() {
		// FIXME it works only for covariates and not for class effect such as species for instance
		matricesZ.clear();
		Matrix matrixZ;
		if (!getRandomEffectsSpecifications().isEmpty()) {
			for (String level : getRandomEffectsSpecifications().keySet()) {
//				int indexOfThisLevel = dataSet.getIndexOfThisField(level);
				List<String> effects = getRandomEffectsSpecifications().get(level);
				matrixZ = new Matrix(getNumberOfObservations(), effects.size());
				for (int i = 0; i < getNumberOfObservations(); i++) {
					for (int j = 0; j < effects.size(); j++) {
						if (effects.get(j).equals("1")) {
							matrixZ.setValueAt(i, j, 1d);
						} else {
							String effectName = effects.get(j);
							int indexOfEffectName = dataSet.getIndexOfThisField(effectName);
							matrixZ.setValueAt(i, j, (Double) dataSet.getValueAt(i, indexOfEffectName));
						}
					}
				}
				matricesZ.put(level, matrixZ);
			}
		}
	}

	
	@Override
	public Set<String> getHierarchicalStructureLevel() {return hierarchicalStructure.keySet();}

	@Override
	public Map<String, DataBlock> getHierarchicalStructure() {return hierarchicalStructure;}

	@Override
	public void setHierarchicalStructureLevel(List<String> hierarchicalStructureLevels) {
		hierarchicalStructure.clear();
		for (int i = 0; i < getNumberOfObservations(); i++) {
			DataBlock currentBlock = null;
			DataBlock parentBlock = null;
			for (int levelIndex = 0; levelIndex < hierarchicalStructureLevels.size(); levelIndex++) {
				String level = hierarchicalStructureLevels.get(levelIndex);
				
				int index = dataSet.getIndexOfThisField(level);	// first check if the field exists
				if (index < 0) {
					throw new InvalidParameterException("Error : This field is not part of the data set : " + level);
				}
				
				String levelFieldValue = dataSet.getValueAt(i, index).toString();

				if (levelIndex == 0) {
					if (!hierarchicalStructure.containsKey(levelFieldValue)) {
						hierarchicalStructure.put(levelFieldValue, new DataBlock(level, levelFieldValue));
					}
					currentBlock = hierarchicalStructure.get(levelFieldValue);
				} else {
					parentBlock = currentBlock;
					if (!parentBlock.containsKey(levelFieldValue)) {
						parentBlock.put(levelFieldValue, new DataBlock(level, levelFieldValue));
					}
					currentBlock = parentBlock.get(levelFieldValue);
				}
				
				currentBlock.addIndex(i);
				
			}
			
		}
	}
	
	
	@Override
	public boolean isThereAnyHierarchicalStructure() {return !hierarchicalStructure.isEmpty();}

	
	@Override
	public Map<String, Matrix> getMatrixZ() {
		return matricesZ;
	}

}
