/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2021 Mathieu Fortin for Rouge-Epicea
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
package repicea.simulation.thinners;

/**
 * The REpiceaTreatmentDefinition interface ensures a particular description of treatment contains 
 * a treatment type and a delay before re-entry.
 * @author Mathieu Fortin - March 2021
 */
public interface REpiceaTreatmentDefinition {

	/**
	 * Provide the treatment type associated with this instance of REpiceaTreatmentDescription.
	 * @return an Enum that stands for the treatment type
	 */
	public Enum getTreatmentType();
	
	/**
	 * Provide the number of years before re-entry. 
	 * @return an integer that is the number of years.
	 */
	public int getDelayBeforeReentryYrs();
	
}
