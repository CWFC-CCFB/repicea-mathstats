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

import java.util.List;
import java.util.Map;
import java.util.Set;

import repicea.math.Matrix;

public interface HierarchicalStatisticalDataStructure extends StatisticalDataStructure {

	/**
	 * This method returns the list of the different levels of the hierarchical structure.
	 * @return a Set of String instances
	 */
	public Set<String> getHierarchicalStructureLevel();
	
	/**
	 * This method returns the index of the observations by levels and id within the level.
	 * @return a Map instance that contains the outer subject levels as keys and DataBlock instances as values
	 */
	public Map<String, DataBlock> getHierarchicalStructure();

	/**
	 * This method set the hierarchical structures of the data set.
	 * @param hierarchicalStructureLevels a List of Strings that contains the field names that serve as hierarchical index
	 * @throws StatisticalDataException
	 */
	public void setHierarchicalStructureLevel(List<String> hierarchicalStructureLevels) throws StatisticalDataException;
	
	/**
	 * Return true if the hierarchical structure has been set or false otherwise.
	 * @return a boolean instance
	 */
	public boolean isThereAnyHierarchicalStructure();

	/**
	 * This method returns a Map of Matrix instance. The key refers to the hierarchical level.
	 * @return a Map instance
	 */
	public Map<String, Matrix> getMatrixZ();
}
