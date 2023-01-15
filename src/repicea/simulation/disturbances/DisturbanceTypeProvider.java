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
package repicea.simulation.disturbances;

import java.util.ArrayList;
import java.util.List;

import repicea.util.REpiceaTranslator;
import repicea.util.REpiceaTranslator.TextableEnum;

public interface DisturbanceTypeProvider {

	public static enum DisturbanceType implements TextableEnum {
		Windstorm("Windstorm","Temp\u00EAte", true),
		Drought("Drought", "S\u00E9cheresse", true),
		Harvest("Harvest", "Coupe", false),
		SpruceBudwormOutbreak("Spruce budworm outbreak", "Epid\u00E9mie de tordeuse des bourgeons de l'\u00E9pinette", true);

		private static List<DisturbanceType> LargeScaleDisturbances;
		
		final boolean largeScale;
		
		DisturbanceType(String englishText, String frenchText, boolean largeScale) {
			setText(englishText, frenchText);
			this.largeScale = largeScale;
		}
		
		
		@Override
		public void setText(String englishText, String frenchText) {
			REpiceaTranslator.setString(this, englishText, frenchText);
		}
		
		@Override
		public String toString() {return REpiceaTranslator.getString(this);}

		/**
		 * This method returns the list of large scale disturbances, mainly for a
		 * proper disturbance occurrence recording implementation. The list is a
		 * copy of an internal list to ensure that the internal list cannot be 
		 * modified.
		 * @return a List of DisturbanceType enum
		 */
		public static List<DisturbanceType> getLargeScaleDisturbances() {
			if (LargeScaleDisturbances == null) {
				LargeScaleDisturbances = new ArrayList<DisturbanceType>();
				for (DisturbanceType type : DisturbanceType.values()) {
					if (type.largeScale) {
						LargeScaleDisturbances.add(type);
					}
				}
			}
			List<DisturbanceType> outputList = new ArrayList<DisturbanceType>();
			outputList.addAll(LargeScaleDisturbances);
			return outputList;
		}
		
	}
	
	/**
	 * Returns the type of disturbance
	 * @return a DisturbanceType enum
	 */
	public DisturbanceType getDisturbanceType();
}
