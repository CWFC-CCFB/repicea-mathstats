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
package repicea.stats.model.glm.copula;

import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

import repicea.math.AbstractMathematicalFunctionWrapper;
import repicea.math.Matrix;
import repicea.stats.LinearStatisticalExpression;
import repicea.stats.data.HierarchicalStatisticalDataStructure;
import repicea.stats.data.StatisticalDataException;
import repicea.stats.model.StatisticalModel;


/**
 * This class is a basic class for all the copula expressions contained in the CopulaLibrary class.
 * @author Mathieu Fortin - October 2011
 */
@SuppressWarnings("serial")
public abstract class CopulaExpression extends AbstractMathematicalFunctionWrapper {
	
	protected List<String> levels = new ArrayList<String>();
	
	protected CopulaExpression(String hierarchicalLevels) {
		super(new LinearStatisticalExpression());
		StringTokenizer tkz = new StringTokenizer(hierarchicalLevels, "/");
		for (int i = 0; i < tkz.countTokens(); i++) {
			levels.add(tkz.nextToken().trim());
		}
	}

	@Override
	public LinearStatisticalExpression getOriginalFunction() {return (LinearStatisticalExpression) super.getOriginalFunction();}
	
	protected List<String> getHierarchicalLevelSpecifications() {return levels;}
	
	protected abstract void setX(int indexFirstObservation, int indexSecondObservation);
	
	protected void initialize(StatisticalModel model, HierarchicalStatisticalDataStructure data) throws StatisticalDataException {
		data.setHierarchicalStructureLevel(levels);
	}
	
}
