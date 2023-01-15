/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2016 Mathieu Fortin for Rouge-Epicea
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
package repicea.simulation;

public class REpiceaPredictorEvent {

	public static class ModelBasedSimulatorEventProperty {
		
		public static final ModelBasedSimulatorEventProperty DEFAULT_BETA_JUST_SET = new ModelBasedSimulatorEventProperty("DEFAULT_BETA_JUST_SET");
		public static final ModelBasedSimulatorEventProperty DEFAULT_RANDOM_EFFECT_AT_THIS_LEVEL_JUST_SET = new ModelBasedSimulatorEventProperty("DEFAULT_RANDOM_EFFECT_AT_THIS_LEVEL_JUST_SET");
		public static final ModelBasedSimulatorEventProperty DEFAULT_RESIDUAL_ERROR_JUST_SET = new ModelBasedSimulatorEventProperty("DEFAULT_RESIDUAL_ERROR_JUST_SET");
		public static final ModelBasedSimulatorEventProperty BLUPS_JUST_SET = new ModelBasedSimulatorEventProperty("BLUPS_JUST_SET");
		public static final ModelBasedSimulatorEventProperty PARAMETERS_DEVIATE_JUST_GENERATED = new ModelBasedSimulatorEventProperty("PARAMETERS_DEVIATE_JUST_GENERATED");
		public static final ModelBasedSimulatorEventProperty RANDOM_EFFECT_DEVIATE_JUST_GENERATED = new ModelBasedSimulatorEventProperty("RANDOM_EFFECT_DEVIATE_JUST_GENERATED");
		public static final ModelBasedSimulatorEventProperty RESIDUAL_ERROR_DEVIATE_JUST_GENERATED = new ModelBasedSimulatorEventProperty("RESIDUAL_ERROR_DEVIATE_JUST_GENERATED");

		private String propertyName;

		protected ModelBasedSimulatorEventProperty(String propertyName) {
			this.propertyName = propertyName;
		}
		
		public String getPropertyName() {return propertyName;}
		
		@Override
		public String toString() {return getPropertyName();}
	}
	
	private final String propertyName;
	private final Object oldValue;
	private final Object newValue;
	private final REpiceaPredictor source;
	
	protected REpiceaPredictorEvent(ModelBasedSimulatorEventProperty property, Object oldValue, Object newValue, REpiceaPredictor source) {
		this.propertyName = property.propertyName;
		this.oldValue = oldValue;
		this.newValue = newValue;
		this.source = source;
	}
		
	
	public String getPropertyName() {return propertyName;}
	public Object getOldValue() {return oldValue;}
	public Object getNewValue() {return newValue;}
	public REpiceaPredictor getSource() {return source;}
	
}
	
	
