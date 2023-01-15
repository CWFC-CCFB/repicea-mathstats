/*
 * This file is part of the repicea-simulation library.
 *
 * Copyright (C) 2009-2014 Mathieu Fortin for Rouge-Epicea
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
package repicea.simulation.processsystem;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Stroke;
import java.util.HashMap;
import java.util.Map;

import javax.swing.BorderFactory;
import javax.swing.ImageIcon;
import javax.swing.border.Border;

import repicea.gui.CommonGuiUtility;
import repicea.simulation.processsystem.ToolPanel.CreateLinkButton;
import repicea.simulation.processsystem.ToolPanel.CreateProcessorButton;
import repicea.simulation.processsystem.ToolPanel.MoveProcessorButton;
import repicea.util.REpiceaTranslator;
import repicea.util.REpiceaTranslator.TextableEnum;

public class UISetup {
	
	protected static enum MessageID implements TextableEnum {
		CreateProcessorButtonToolTip("Create a processor", "Cr\u00E9er un processeur"),
		CreateLinkButtonToolTip("Create a link between two processors", "Cr\u00E9er un lien entre deux processeurs"),
		MoveProcessorButtonToolTip("Move a processor", "D\u00E9placer un processeur");

		MessageID(String englishText, String frenchText) {
			setText(englishText, frenchText);
		}
		
		@Override
		public void setText(String englishText, String frenchText) {
			REpiceaTranslator.setString(this, englishText, frenchText);
		}

		@Override
		public String toString() {return REpiceaTranslator.getString(this);}
	}
	
	
	
	protected static Border ButtonSelectedBorder = BorderFactory.createLineBorder(Color.BLACK, 3);
	
	public static Border ButtonDefaultBorder = BorderFactory.createLineBorder(Color.BLACK, 1);
	
	public static Stroke BoldStroke = new BasicStroke(4);

	public static Stroke DefaultStroke = new BasicStroke(2);
	
	public static int XGap = 100;
	
	public static int YGap = 100;
	
	protected static int XOrigin = 100;
	
	protected static int YOrigin = 100;
	
	protected static enum BasicMode {CreateLink, MoveProcessor};
	
	protected static Color DefaultColor = Color.BLACK;

	public static final Map<String, ImageIcon> Icons = new HashMap<String, ImageIcon>();
	static {
		Icons.put(ProcessorButton.class.getName(), CommonGuiUtility.retrieveIcon(UISetup.class, "sawmill.png"));
		Icons.put("LinkButtonIcon",  CommonGuiUtility.retrieveIcon(UISetup.class, "linkIcon.png"));
		Icons.put("MoveButton", CommonGuiUtility.retrieveIcon(UISetup.class, "handPointerIcon.png"));
	}
	
	public static final Map<String, String> ToolTips = new HashMap<String, String>();
	static {
		ToolTips.put(CreateProcessorButton.class.getName(), MessageID.CreateProcessorButtonToolTip.toString());
		ToolTips.put(CreateLinkButton.class.getName(), MessageID.CreateLinkButtonToolTip.toString());
		ToolTips.put(MoveProcessorButton.class.getName(), MessageID.MoveProcessorButtonToolTip.toString());
	}
	

	
//	private static ImageIcon retrieveIcon(String iconName) {
//		String iconPath = ObjectUtility.getRelativePackagePath(UISetup.class) + iconName;
//		InputStream iconInputStream = ClassLoader.getSystemResourceAsStream(iconPath);
//		try {
//			Image image = ImageIO.read(iconInputStream);
//			return new ImageIcon(image);
//		} catch (Exception e2) {
//			return null;
//		}
//	}
	

	
	
}
