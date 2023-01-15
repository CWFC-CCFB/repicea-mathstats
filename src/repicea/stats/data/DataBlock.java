/*
 * This file is part of the repicea library.
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

@SuppressWarnings("serial")
public class DataBlock extends LinkedHashMap<String, DataBlock> {
	
	private String hierarchicalLevel;
	private String subject;
	private DataBlock parentBlock;
	private List<Integer> indices;
	
	protected DataBlock(String hierarchicalLevel, String subject) {
		this.hierarchicalLevel = hierarchicalLevel;
		this.subject = subject;
		indices = new ArrayList<Integer>();
	}
	
	void addIndex(int i) {
		indices.add(i);
	}

	/**
	 * This method returns the list of block corresponding to a particular hierarchical
	 * level.
	 * @param hierarchicalLevel a String that defines the hierarchical level
	 * @return a List of DataBlock instance, which can be empty if no block of this hierarchical level has been found
	 */
	public List<DataBlock> getBlocksOfThisLevel(String hierarchicalLevel) {
		List<DataBlock> outputList = new ArrayList<DataBlock>();
		if (getHierarchicalLevel().equals(hierarchicalLevel)) {
			outputList.add(this);
		} else {
			for (DataBlock db : this.values()) {
				outputList.addAll(db.getBlocksOfThisLevel(hierarchicalLevel));
			}
		}
		return outputList;
	}

	
	private void setParentBlock(DataBlock db) {
		if (parentBlock != null) {
			throw new InvalidParameterException("DataBlock.setParentParent : block has already been set!");
		} else {
			parentBlock = db;
		}
	}
	
	public List<Integer> getIndices() {return indices;}

	@Override
	public DataBlock put(String subjectID, DataBlock db) {
		db.setParentBlock(this);
		return super.put(subjectID, db);
	}
	
	@Override
	public String toString() {
		if (isEmpty()) {
			return hierarchicalLevel + " - " + subject + " " + indices.toString();
		} else {
			return hierarchicalLevel + " - " + subject + " " + super.toString();
		}
	}
	
	


	/**
	 * This method returns the hierarchical level of the data block.
	 * @return a String
	 */
	public String getHierarchicalLevel() {
		return hierarchicalLevel;
	}
	
}
