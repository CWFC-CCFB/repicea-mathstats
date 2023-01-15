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

import java.util.List;

import repicea.simulation.REpiceaBinaryEventPredictor;

@SuppressWarnings("serial")
public abstract class REpiceaThinner<S, T> extends REpiceaBinaryEventPredictor<S, T> {

	protected REpiceaThinner(boolean isParametersVariabilityEnabled, boolean isRandomEffectsVariabilityEnabled,	boolean isResidualVariabilityEnabled) {
		super(isParametersVariabilityEnabled, isRandomEffectsVariabilityEnabled, isResidualVariabilityEnabled);
	}

	/**
	 * Return some information about the treatment.
	 * @param stand a S-derived instance
	 * @return an REpiceaTreatmentDefinition or null if no information can be provided.
	 */
	public abstract REpiceaTreatmentDefinition getTreatmentDefinitionForThisHarvestedStand(S stand);
	
	/**
	 * Return the list of all possible treatments.
	 * @return a List of Enum variables
	 */
	public abstract List<Enum> getTreatmentList();
}
