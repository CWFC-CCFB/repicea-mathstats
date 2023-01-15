/*
 * This file is part of the repicea-foresttools library.
 *
 * Copyright (C) 2009-2013 Mathieu Fortin for Rouge-Epicea
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed with the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * Please see the license at http://www.gnu.org/copyleft/lesser.html.
 */
package repicea.simulation;

import java.io.IOException;

import repicea.io.javacsv.CSVReader;
import repicea.math.Matrix;

/**
 * This class offers static methods to load parameters. The methods rely on a CSVReader instance
 * to read parameter in .csv files 
 * @author Mathieu Fortin - July 2012
 */
public class ParameterLoader {
	
	/**
	 * This method reads a file and retrieve a vector of parameters.
	 * @param filename the path of the file to be read (*.csv)
	 * @return a ParameterMap instance
	 * @throws IOException if something goes wrong while reading the file.
	 */
	public static ParameterMap loadVectorFromFile(int numberOfIndices, String filename) throws IOException {
		ParameterMap pm = new ParameterMap();
//		List<Double> beta = new ArrayList<Double>();
		CSVReader reader = null;
		int index1;
		int index2;
		double parameter;
		try {			
			reader = new CSVReader(filename);
			Object[] lineRead = reader.nextRecord();
			while (lineRead != null) {
				if (numberOfIndices == 2) {
					index1 = Integer.parseInt(lineRead[0].toString());
					index2 = Integer.parseInt(lineRead[1].toString());
					parameter = Double.parseDouble(lineRead[2].toString());
					pm.addParameter(index1, index2, parameter);
				} else if (numberOfIndices == 1) {
					index2 = Integer.parseInt(lineRead[0].toString());
					parameter = Double.parseDouble(lineRead[1].toString());
					pm.addParameter(index2, parameter);
				} else {
					parameter = Double.parseDouble(lineRead[0].toString());
					pm.addParameter(parameter);
				}
				lineRead = reader.nextRecord();
			}
			reader.close();
			return pm;
		} catch (Exception e) {
			if (reader != null) {
				reader.close();
			}
			throw new IOException("ParameterLoader.loadVectorFormFile() : Unable to read table");
		}
	}

	
	/**
	 * This method reads a file and retrieve a vector of parameters.
	 * @param filename the path of the file to be read (*.csv)
	 * @return a ParameterMap instance
	 * @throws IOException if something goes wrong while reading the file.
	 */
	public static ParameterMap loadVectorFromFile(String filename) throws IOException {
		return loadVectorFromFile(0, filename);
	}
	
	
	/**
	 * This method reads the file and retrieve a matrix of parameters.
	 * @param filename the path of the file to be read (*.csv)
	 * @return a Matrix instance
	 * @throws IOException if something goes wrong while reading the file
	 */
	public static Matrix loadMatrixFromFile(String filename) throws IOException {
		CSVReader reader = null;
		try {			
			reader = new CSVReader(filename);
			Matrix omega = new Matrix(reader.getRecordCount(), reader.getFieldCount());
			Object[] lineRead = reader.nextRecord();
			int record = 0;
			while (lineRead != null) {
				for (int i = 0; i < reader.getFieldCount(); i++) {
					double parameter = Double.parseDouble(lineRead[i].toString());
					omega.setValueAt(record, i, parameter);
				}
				lineRead = reader.nextRecord();
				record++;
			}
			reader.close();
			return omega;
		} catch (Exception e) {
			if (reader != null) {
				reader.close();
			}
			throw new IOException("ParameterLoader.loadMatrixFromFile() : Unable to read table");
		}
	}

}
