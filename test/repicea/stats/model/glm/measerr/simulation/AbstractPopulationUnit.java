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
package repicea.stats.model.glm.measerr.simulation;

import java.util.ArrayList;
import java.util.List;

abstract class AbstractPopulationUnit implements PopulationUnit {

	protected static final List<String> FIELDNAMES_INPUT_DATASET = new ArrayList<String>();
	static {
		FIELDNAMES_INPUT_DATASET.add("y");
		FIELDNAMES_INPUT_DATASET.add("distanceToConspecific");
		FIELDNAMES_INPUT_DATASET.add("variance");
		FIELDNAMES_INPUT_DATASET.add("trueDistanceToConspecific");
	}

	
	
	final int id;
	final boolean isConspecificIn;
	double trueDistanceToConspecific;
	double measuredDistanceToConspecific;
	private int y;

	AbstractPopulationUnit(int id, boolean isConspecificIn) {
		this.id = id;
		this.isConspecificIn = isConspecificIn;
	}
	
	@Override
	public abstract Object[] asObservation(boolean detailed);

	@Override
	public abstract PopulationUnit clone();
	
	protected void setY(int y) {this.y = y;}
	
	protected int getY() {return this.y;}
}
