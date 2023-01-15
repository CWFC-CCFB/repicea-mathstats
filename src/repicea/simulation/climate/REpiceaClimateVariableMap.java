/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2019 Mathieu Fortin for Rouge-Epicea
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
package repicea.simulation.climate;

import java.util.HashMap;

import repicea.simulation.climate.REpiceaClimateVariableMap.ClimateVariable;
import repicea.util.REpiceaTranslator;
import repicea.util.REpiceaTranslator.TextableEnum;

@SuppressWarnings("serial")
public class REpiceaClimateVariableMap extends HashMap<ClimateVariable, Double> {
	
	public static enum ClimateVariable implements TextableEnum {
		/**
		 * Mean annual temperature (C)
		 */
		MeanAnnualTempC("Mean annual temperature (\u00B0C)", "Temp\u00E9rature annuelle moyenne (\u00B0C)"),
		
		/**
		 * Mean annual precipitation (mm)
		 */
		MeanAnnualPrecMm("Mean annual precipitation (mm)", "Pr\u00Ecipitations annuelles moyennes (mm)"),
		
		/**
		 * Mean temperature of the growing season (C)
		 */
		MeanGrowingSeasonTempC("Mean growing season temperature (\u00B0C)", "Moyenne des temp\u00E9rature de la saison de croissance (\u00B0C)"),
		
		/**
		 * Mean precipitation of the growing season (mm)
		 */
		MeanGrowingSeasonPrecMm("Mean growing season precipitation (mm)", "Moyenne des pr\u00E9cipitations de la saison de croissance (mm)");

		ClimateVariable(String englishText, String frenchText) {
			setText(englishText, frenchText);
		}
		
		@Override
		public void setText(String englishText, String frenchText) {
			REpiceaTranslator.setString(this, englishText, frenchText);
		}
		
		@Override
		public String toString() {return REpiceaTranslator.getString(this);}
	}


	/**
	 * This interface ensures that the REpiceaClimateVariableMap instance
	 * can be updated according to climate change forecasts. 
	 * @author Mathieu Fortin - June 2019
	 */
	public static interface UpdatableClimateVariableMap {
		
		/**
		 * Performs the change on the climate variable following a linear extrapolation pattern.
		 * @param climateChangeTrend a REpiceaClimateChangeTrend instance
		 * @param dateYr the current date (years)
		 */
		public REpiceaClimateVariableMap getUpdatedClimateVariableMap(REpiceaClimateChangeTrend climateChangeTrend, int dateYr);
		
	}
	
}
