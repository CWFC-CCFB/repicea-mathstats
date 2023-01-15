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

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;

import javax.swing.BorderFactory;
import javax.swing.Icon;
import javax.swing.ImageIcon;

import repicea.gui.permissions.REpiceaGUIPermission;

@SuppressWarnings("serial")
public abstract class SelectableJButton extends AbstractPermissionProviderButton {

	protected final ImageIcon originalIcon;
	private int borderWidth = 1;
	private Color borderColor = Color.BLACK;
	
	protected SelectableJButton(REpiceaGUIPermission permission) {
		super(permission);
		
		String className = getClass().getName();
		Icon icon = UISetup.Icons.get(className);
		if (icon == null) {
			icon = getDefaultIcon();	// default value
		}
		setIcon(icon);
		originalIcon = (ImageIcon) icon;
		setToolTip();
		Dimension dim;
		if (getIcon() != null) {
			dim = new Dimension(originalIcon.getIconWidth(), originalIcon.getIconHeight());
		} else {
			dim = new Dimension(30,30);
		}
		setPreferredSize(dim);
		setMinimumSize(dim);
		setMaximumSize(dim);
		setSize(dim);
		setBorder(UISetup.ButtonDefaultBorder);
		setFocusable(false);
	}

	private void setToolTip() {
		String toolTip = UISetup.ToolTips.get(getClass().getName());
		if (toolTip != null) {
			setToolTipText(toolTip);
		}
	}

	protected abstract Icon getDefaultIcon();
	
	@Override
	public void setSelected(boolean bool) {
		super.setSelected(bool);
		if (bool) {
			setBorderWidth(3);
		} else {
			setBorderWidth(1);
		}
	}

//	@Deprecated
//	public void setBorder(Border border) {
//		super.setBorder(border);
//	}
	
	
	protected void setBorderWidth(int width) {
		if (width != borderWidth) {
			borderWidth = width;
			setBorder(BorderFactory.createLineBorder(borderColor, borderWidth));
		}
	}
	
	protected void setBorderColor(Color col) {
		if (!col.equals(borderColor)) {
			borderColor = col;
			setBorder(BorderFactory.createLineBorder(borderColor, borderWidth));
		}
	}
		
	@Override
	public void paint(Graphics g) {
//		setBorder(BorderFactory.createLineBorder(borderColor, borderWidth));
		super.paint(g);
	}
}
