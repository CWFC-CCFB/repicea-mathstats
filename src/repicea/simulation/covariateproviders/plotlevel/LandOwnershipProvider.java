/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2017 Mathieu Fortin for Rouge-Epicea
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
package repicea.simulation.covariateproviders.plotlevel;

import java.util.HashMap;
import java.util.Map;

public interface LandOwnershipProvider {
	
	public static enum LandOwnership {
		Public,
		Private;
		
		private static Map<String, LandOwnership> OwnershipMap;
		
		/**
		 * This method returns the LandOwnership associated with the ownership code, which can be "PU" for public or "PR"
		 * for private-owned land.
		 * @param ownershipCode a String
		 * @return a LandOwnership enum
		 */
		public static LandOwnership getLandOwnership(String ownershipCode) {
			if (OwnershipMap == null) {
				OwnershipMap = new HashMap<String, LandOwnership>();
				for (LandOwnership landOwnership : LandOwnership.values()) {
					OwnershipMap.put(landOwnership.name().substring(0, 2).toUpperCase(), landOwnership);
				}
			}
			return (OwnershipMap.get(ownershipCode));
		}
		
		
		
	}

	
	/**
	 * This method returns the land ownership of the plot instance.
	 * @return a LandOwnership instance
	 */
	public LandOwnership getLandOwnership();
	
	
}
