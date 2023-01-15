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

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.util.ArrayList;
import java.util.List;

import javax.swing.AbstractButton;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.Icon;
import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.JSeparator;
import javax.swing.KeyStroke;
import javax.swing.SwingConstants;

import repicea.gui.REpiceaPanel;
import repicea.gui.Resettable;
import repicea.gui.permissions.REpiceaGUIPermission;
import repicea.gui.permissions.REpiceaGUIPermissionProvider;
import repicea.simulation.processsystem.UISetup.BasicMode;

@SuppressWarnings("serial")
public class ToolPanel extends REpiceaPanel implements ActionListener, Resettable, REpiceaGUIPermissionProvider {

	protected class KeyEventHandler implements ActionListener {

		protected final ToolButton button;
		
		protected KeyEventHandler(ToolButton button) {
			this.button = button;
		}
		
		@Override
		public void actionPerformed(ActionEvent arg0) {
			button.doClick();
		}
		
	}
	
	protected class CreateProcessorButton extends DnDCompatibleButton {

		protected CreateProcessorButton(REpiceaGUIPermission permission) {
			super(permission);
		}
		
		@Override
		public Processor createNewProcessor() {
			return new Processor();
		}

		@Override
		protected Icon getDefaultIcon() {
			return UISetup.Icons.get(ProcessorButton.class.getName());
		}
	}
	

	protected class CreateLinkButton extends ToolButton {

		protected CreateLinkButton(REpiceaGUIPermission permission) {
			super(permission);
			mode = BasicMode.CreateLink;
		}
		
		@Override
		protected Icon getDefaultIcon() {
			return UISetup.Icons.get("LinkButtonIcon");
		}

	}

	protected class MoveProcessorButton extends ToolButton {

		protected MoveProcessorButton(REpiceaGUIPermission permission) {
			super(permission);
			mode = BasicMode.MoveProcessor;
		}

		@Override
		protected Icon getDefaultIcon() {
			return UISetup.Icons.get("MoveButton");
		}

	}

	private final static List<Integer> keyStrokes = new ArrayList<Integer>();
	static {
		keyStrokes.add(KeyEvent.VK_F1);
		keyStrokes.add(KeyEvent.VK_F2);
		keyStrokes.add(KeyEvent.VK_F3);
		keyStrokes.add(KeyEvent.VK_F4);
		keyStrokes.add(KeyEvent.VK_F5);
		keyStrokes.add(KeyEvent.VK_F6);
		keyStrokes.add(KeyEvent.VK_F7);
		keyStrokes.add(KeyEvent.VK_F8);
		keyStrokes.add(KeyEvent.VK_F9);
		keyStrokes.add(KeyEvent.VK_F10);
	}
	
	private final List<ToolButton> selectableButtons;
	protected final List<DnDCompatibleButton> dndButtons;

	private DnDCompatibleButton createProcessorButton;
	private ToolButton createLinkButton;
	protected ToolButton simpleSelectionButton;

	
	protected final SystemPanel owner;
	
	protected ToolPanel(SystemPanel owner) {
		super();
		this.owner = owner;

		dndButtons = new ArrayList<DnDCompatibleButton>();
		selectableButtons = new ArrayList<ToolButton>();
		init();
		createUI();
	}

	
	protected void init() {
		createProcessorButton = new CreateProcessorButton(getGUIPermission());
		dndButtons.add(createProcessorButton);

		simpleSelectionButton = new MoveProcessorButton(getGUIPermission());
		simpleSelectionButton.setName("simpleSelectionButton");		
		addSelectableButton(simpleSelectionButton);
		simpleSelectionButton.setSelected(true);
		
		createLinkButton = new CreateLinkButton(getGUIPermission());
		createLinkButton.setName("createLinkButton");
		addSelectableButton(createLinkButton);
	}
	
	protected void addSelectableButton(ToolButton button) {
		selectableButtons.add(button);
		if (selectableButtons.size() <= keyStrokes.size()) {
			KeyStroke stroke = KeyStroke.getKeyStroke(keyStrokes.get(selectableButtons.size() - 1), 0);
			registerKeyboardAction(new KeyEventHandler(button), stroke, JComponent.WHEN_IN_FOCUSED_WINDOW);
		}
	}
	
	
	private JPanel putInAJPanel(Component comp) {
		JPanel pane = new JPanel(new FlowLayout(FlowLayout.CENTER));
		pane.add(comp);
		return pane;
	}
	
	protected void createUI() {
		setLayout(new BorderLayout());
		JPanel panel = new JPanel();
		panel.setLayout(new BoxLayout(panel, BoxLayout.Y_AXIS));
		
		panel.add(Box.createVerticalStrut(5));
		for (AbstractButton button : dndButtons) {
			panel.add(putInAJPanel(button));
			panel.add(Box.createVerticalStrut(5));
		}
		panel.add(new JSeparator(SwingConstants.HORIZONTAL));
		
		panel.add(Box.createVerticalStrut(5));
		for (SelectableJButton button : selectableButtons) {
			panel.add(putInAJPanel(button));
			panel.add(Box.createVerticalStrut(5));
		}
		
		panel.add(Box.createVerticalGlue());
		
		add(panel, BorderLayout.NORTH);
		setBorder(UISetup.ButtonDefaultBorder);
	}
	
	@Override
	public void refreshInterface() {
		for (ToolButton b : selectableButtons) {
			if (b.getMode().equals(owner.getMode())) {
				b.setSelected(true);
				break;
			}
		}
	}

	@Override
	public void listenTo() {
		for (SelectableJButton button : selectableButtons) {
			button.addActionListener(this);
			button.addActionListener(owner);
		}
	}

	@Override
	public void doNotListenToAnymore() {
		for (SelectableJButton button : selectableButtons) {
			button.removeActionListener(this);
			button.removeActionListener(owner);
		}
	}



	@Override
	public void actionPerformed(ActionEvent arg0) {
		if (arg0.getSource() instanceof SelectableJButton) {
			SelectableJButton selectedButton = (SelectableJButton) arg0.getSource();
			setSelectedButton(selectedButton);
		}
	}

	private void setSelectedButton(SelectableJButton selectedButton) {
		for (SelectableJButton button : selectableButtons) {
			button.setSelected(selectedButton.equals(button));
		}
	}

	ToolButton getSelectedButton() {
		for (ToolButton button : selectableButtons) {
			if (button.isSelected()) {
				return button;
			}
		}
		return null;
	}

	@Override
	public void reset() {
		doNotListenToAnymore();
		for (SelectableJButton button : selectableButtons) {
			button.setSelected(false);
		}
		listenTo();
	}


	@Override
	public REpiceaGUIPermission getGUIPermission() {
		return owner.getListManager().getGUIPermission();
	}

//	@Override
//	public void keyPressed(KeyEvent arg0) {
//		if (keyStrokes.contains(arg0.getKeyCode())) {
//			int index = keyStrokes.indexOf(arg0.getKeyCode());
//			if (index < selectableButtons.size()) {
//				selectableButtons.get(index).setSelected(true);
//			}
//		}
//	}
//
//
//	@Override
//	public void keyReleased(KeyEvent arg0) {
//		int u = 0;
//	}
//
//
//	@Override
//	public void keyTyped(KeyEvent arg0) {
//		int u = 0;
//	}


}
