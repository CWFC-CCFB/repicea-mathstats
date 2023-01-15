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

import repicea.math.Matrix;
import repicea.stats.StatisticalUtility;

class SpatialPopulationUnit extends AbstractPopulationUnit {

	private static final List<String> FIELDNAMES = new ArrayList<String>();
	static {
		FIELDNAMES.add("y");
		FIELDNAMES.add("xCoord");
		FIELDNAMES.add("yCoord");
		FIELDNAMES.add("distanceToConsp");
		FIELDNAMES.add("isConspecificIn");
	}

	
	
	final int xCoord;
	final int yCoord;
	final DistanceCalculator dc;
	double variance;

	SpatialPopulationUnit(int id, int xCoord, int yCoord, boolean isConspecificIn) {
		super(id, isConspecificIn);
		this.xCoord = xCoord;
		this.yCoord = yCoord;
		this.dc = new DistanceCalculator(this);
	}
	
	double getDistanceFrom(SpatialPopulationUnit thatUnit) {
		double xDiff = xCoord - thatUnit.xCoord;
		double yDiff = yCoord - thatUnit.yCoord;
		return Math.sqrt(xDiff * xDiff + yDiff * yDiff);
	}
	
	@Override
	public int setY(Matrix parameters) {
		double pred = 1 - Math.exp(-Math.exp(parameters.getValueAt(0, 0) + parameters.getValueAt(1, 0) * trueDistanceToConspecific));
		boolean occurred = StatisticalUtility.getRandom().nextDouble() < pred;
		setY(occurred ? 1 : 0);
		return getY();
	}

	@Override
	public Object[] asObservation(boolean detailed) {
		if (detailed) {
			Object[] record = new Object[5];
			record[0] = getY();
			record[1] = xCoord;
			record[2] = yCoord;
			record[3] = trueDistanceToConspecific;
			record[4] = isConspecificIn;
			return record;
		} else {
			Object[] record = new Object[4];
			record[0] = getY();
			record[1] = measuredDistanceToConspecific;
			record[2] = variance;
			record[3] = trueDistanceToConspecific;
			return record;
		}
	}

	@Override
	public SpatialPopulationUnit clone() {
		SpatialPopulationUnit pu = new SpatialPopulationUnit(this.id, this.xCoord, this.yCoord, this.isConspecificIn);
		pu.trueDistanceToConspecific = trueDistanceToConspecific;
		pu.measuredDistanceToConspecific = measuredDistanceToConspecific;
		pu.setY(getY());
		return pu;
	}

	@Override
	public List<String> getFieldname(boolean detailed) {
		if (detailed) 
			return FIELDNAMES;
		else 
			return AbstractPopulationUnit.FIELDNAMES_INPUT_DATASET;
	}

}
