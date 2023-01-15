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

import java.security.InvalidParameterException;

import repicea.simulation.climate.REpiceaClimateGenerator.ClimateChangeOption;
import repicea.simulation.climate.REpiceaClimateGenerator.ClimateChangeScenario;
import repicea.simulation.climate.REpiceaClimateGenerator.RepresentativeConcentrationPathway;
import repicea.simulation.covariateproviders.plotlevel.GeographicalCoordinatesProvider;
import repicea.util.REpiceaTranslator;
import repicea.util.REpiceaTranslator.TextableEnum;

/**
 * The REpiceaClimateGenerator interface ensures that the climate generator returns 
 * a REpiceaClimateVariableMap instance.
 * @author Mathieu Fortin - June 2019
 */
public interface REpiceaClimateGenerator<P extends GeographicalCoordinatesProvider> {

	/**
	 * An empty interface used to identify climate change options.
	 * @author Mathieu Fortin - September 2021
	 */
	public interface ClimateChangeScenario extends TextableEnum {}
	
	public interface ClimateChangeScenarioProvider {
		ClimateChangeScenario getClimageChangeScenario();
	}
	
	
	public static enum RepresentativeConcentrationPathway implements ClimateChangeScenario {
		RCP2_6("RCP 2.6", "RCP 2.6"),
		RCP4_5("RCP 4.5", "RCP 4.5"),
		RCP6_0("RCP 6.0", "RCP 6.0"),
		RCP8_5("RCP 8.5", "RCP 8.5");

		RepresentativeConcentrationPathway(String englishText, String frenchText) {
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
	 * This enum serves to set the climate change options in the initial parameters
	 * @author Mathieu Fortin - December 2010
	 */
	@Deprecated	
	public static enum ClimateChangeOption implements ClimateChangeScenario {
		NoChange("No change", 
				"Aucun changement"), 
		Plus2Degrees("2-degree increase over the XXI century", 
				"Augmentation de 2 degr\u00E9s au cours du 21e si\u00E8cle"), 
		Plus4Degrees("4-degree increase over the XXI century", 
				"Augmentation de 4 degr\u00E9s au cours du 21e si\u00E8cle"), 
		Plus6Degrees("6-degree increase over the XXI century", 
				"Augmentation de 6 degr\u00E9s au cours du 21e si\u00E8cle");

		ClimateChangeOption(String englishText, String frenchText) {
			setText(englishText, frenchText);
		}
		
		@Override
		public void setText(String englishText, String frenchText) {
			REpiceaTranslator.setString(this, englishText, frenchText);
		}
		
		@Override
		public String toString() {
			return REpiceaTranslator.getString(this);
		}		
	}

	/**
	 * Returns a map of climate variables depending on the geographical coordinates of the plot.
	 * @param plot a GeographicalCoordinatesProvider instance
	 * @return a REpiceaClimateVariableMap-derived instance
	 */
	public REpiceaClimateVariableMap getClimateVariables(P plot);
	
	public static class ClimateChangeScenarioHelper {
		@SuppressWarnings("deprecation")
		static public ClimateChangeScenario getClimateChangeScenarioFromString(String climateChangeScenarioStr) {
			try {
				return RepresentativeConcentrationPathway.valueOf(climateChangeScenarioStr);
			} catch(Exception e) {
				try {
					return ClimateChangeOption.valueOf(climateChangeScenarioStr);
				} catch (Exception e2) {
					throw new InvalidParameterException("The climate change scenario " + climateChangeScenarioStr + " is not recognized!");
				}
			}
		}
	}
}
