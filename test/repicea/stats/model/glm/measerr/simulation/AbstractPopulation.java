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

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import repicea.io.FormatField;
import repicea.io.javacsv.CSVField;
import repicea.io.javacsv.CSVWriter;
import repicea.math.Matrix;

abstract class AbstractPopulation<P extends PopulationUnit> implements Population<P> {

	protected static final List<Double> POTENTIAL_X_VALUES = new ArrayList<Double>();
	static {
		for (double d = 0; d < 14.99; d = d + 0.1) {
			POTENTIAL_X_VALUES.add(d);
		}
		for (double d = 15; d < 900; d = d + 1) {
			POTENTIAL_X_VALUES.add(d);
		}
	}

	
	final List<P> populationUnits;
	final Matrix trueBeta;
	final int size;
	
	AbstractPopulation(Matrix trueBeta, int size) {
		this.trueBeta = trueBeta.getDeepClone();
		this.size = size;
		this.populationUnits = new ArrayList<P>();
	}
	
	@Override
	public final int generateRealizations() {
		int pos = 0;
		for (P pu : populationUnits) {
			pos += pu.setY(trueBeta);
		}
		return pos;
	}

	private List<FormatField> getFormatFields(List<String> fieldnames) {
		List<FormatField> fields = new ArrayList<FormatField>();
		for (String s : fieldnames) {
			fields.add(new CSVField(s));
		}
		return fields;
	}
	
	void save(String filename) throws IOException {
		CSVWriter writer = new CSVWriter(new File(filename), false);
		writer.setFields(getFormatFields(populationUnits.get(0).getFieldname(true)));	// true: detailed field names
		
		for (PopulationUnit pu : this.populationUnits) {
			writer.addRecord(pu.asObservation(true));
		}
		writer.close();
	}

	@Override
	public abstract Population<P> clone();
}
