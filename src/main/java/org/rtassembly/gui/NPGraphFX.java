/*****************************************************************************
 * Copyright (c) Son Hoang Nguyen, IMB - UQ, All rights reserved.         *
 *                                                                           *
 * Redistribution and use in source and binary forms, with or without        *
 * modification, are permitted provided that the following conditions        *
 * are met:                                                                  * 
 *                                                                           *
 * 1. Redistributions of source code must retain the above copyright notice, *
 *    this list of conditions and the following disclaimer.                  *
 * 2. Redistributions in binary form must reproduce the above copyright      *
 *    notice, this list of conditions and the following disclaimer in the    *
 *    documentation and/or other materials provided with the distribution.   *
 * 3. Neither the names of the institutions nor the names of the contributors*
 *    may be used to endorse or promote products derived from this software  *
 *    without specific prior written permission.                             *
 *                                                                           *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS   *
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, *
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR    *
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR         *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,     *
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       *
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        *
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      *
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              *
 ****************************************************************************/

/**************************     REVISION HISTORY    **************************
 * 03/25/2018 - Son Hoang Nguyen: Created                                        
 *  
 ****************************************************************************/
package org.rtassembly.gui;

import java.io.File;
import java.io.IOException;
import java.util.concurrent.Executors;
import java.util.concurrent.ScheduledExecutorService;

import org.graphstream.stream.thread.ThreadProxyPipe;
import org.graphstream.ui.fx_viewer.FxDefaultView;
import org.graphstream.ui.fx_viewer.FxViewer;
import org.graphstream.ui.javafx.FxGraphRenderer;
import org.graphstream.ui.view.Viewer;
import org.graphstream.ui.view.Viewer.CloseFramePolicy;
import org.rtassembly.npgraph.Alignment;
import org.rtassembly.npgraph.GraphUtil;
import org.rtassembly.npgraph.HybridAssembler;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import japsa.util.FxDialogs;
import japsa.util.ImageButton;
import javafx.application.Application;
import javafx.application.Platform;
import javafx.geometry.HPos;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.Scene;
import javafx.scene.control.Button;
import javafx.scene.control.CheckBox;
import javafx.scene.control.ComboBox;
import javafx.scene.control.Label;
import javafx.scene.control.Separator;
import javafx.scene.control.Tab;
import javafx.scene.control.TabPane;
import javafx.scene.control.TextField;
import javafx.scene.image.Image;
import javafx.scene.image.ImageView;
import javafx.scene.input.KeyCode;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.ColumnConstraints;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.RowConstraints;
import javafx.scene.layout.StackPane;
import javafx.scene.layout.VBox;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.scene.text.Text;
import javafx.stage.DirectoryChooser;
import javafx.stage.FileChooser;
import javafx.stage.FileChooser.ExtensionFilter;
import javafx.stage.Stage;


@SuppressWarnings("restriction")
public class NPGraphFX extends Application{
    private static final Logger LOG = LoggerFactory.getLogger(NPGraphFX.class);

	static HybridAssembler myass = new HybridAssembler();
	static Viewer graphViewer;
	
	public static void setAssembler(HybridAssembler ass){
		myass = ass;
	}
    
    public void start(Stage primaryStage){  	 
    	Stage graphStage = new Stage();
    	
        try {
        	/*
        	 * Setup graph stage here
        	 */
            BorderPane gBorder = new BorderPane();
            gBorder.setCenter(addGraphResolverPane());
            Scene gscene = new Scene(gBorder);
            graphStage.setScene(gscene);
            graphStage.setTitle("Assembly graph");
//            graphStage.show();
            
            
            /*
             * Setup primary stage here
             */
        	//Start with a BorderPane as root
    	    BorderPane pBorder = new BorderPane();
    	    // Put start/stop button here  
    	    HBox hbox = addHBox();
    	    pBorder.setTop(hbox);
    	    // All the parameters setting to the left 
    	    //TODO: pass graphStage here to bind with graph file being loaded
    	    leftBox=addVBox(primaryStage, graphStage);
    	    pBorder.setLeft(leftBox);
            
            // Here the main content    
            pBorder.setCenter(addAssemblyStatsPane());
            
            Scene pscene = new Scene(pBorder);
            primaryStage.setScene(pscene);
            primaryStage.setTitle("Main");
            primaryStage.setOnCloseRequest(e -> {
                Platform.exit();
                System.exit(0);
            });
            primaryStage.show();
        	
            
            /*
             * A thread to run the algorithm
             */
    	    new Thread(new Runnable(){

    			@Override
    			public void run() {
    				while (!myass.getReady()){
    					//LOG.info("NOT READY");
    					try {
    						Thread.sleep(1000);
    					} catch (InterruptedException e) {					
    						e.printStackTrace();
    					}			
    				}
    				LOG.info("GO");
//    				updateData();
    				try{
    					myass.assembly();
    					myass.postProcessGraph();
    				}catch (Exception e){
    					System.err.println(e.getMessage());
    					e.printStackTrace();
    				}

    			}
    			
    		}).start();
            
        } catch (Exception exc) {
            exc.printStackTrace();
        }
    	
    }
    
    /*
     * Components from left pane
     */
    private VBox leftBox;
    private TabPane tabPane; 
    private Tab statsTab, graphTab;
    private Button 	buttonStart, buttonStop, buttonRestart,
    				shortInputBrowseButton, longReadsBrowseButton, outputBrowseButton, mm2BrowseButton;
    private TextField shortInputTF, longInputTF, outputTF, minQualTF, mm2PathTF, mm2OptTF;
    private CheckBox graphCB, overwriteCB;
    private ComboBox<String> shortInputFormatCombo, longInputFormatCombo;
    /*
     * Creates an HBox with two buttons for the top region
     */
        
    private HBox addHBox() {
 
        HBox hbox = new HBox();
        hbox.setPadding(new Insets(15, 12, 15, 12));
        hbox.setSpacing(10);   // Gap between nodes
        hbox.setStyle("-fx-background-color: #336699;");
 

        Image imageStart = new Image(getClass().getResourceAsStream("/start.png"));
        ImageView viewStart = new ImageView(imageStart); 
        viewStart.setFitWidth(20);
        viewStart.setFitHeight(20);
        buttonStart = new Button("Start", viewStart);
        buttonStart.setPrefSize(100, 20);
        buttonStart.setOnAction((event) -> {
        	if(!checkingAndSetting())
        		return;
        	
			//Start running
			leftBox.setDisable(true);
			//TODO: leftbox slide away
			buttonStart.setDisable(true);;
			buttonStop.setDisable(false);
		});
        
        Image imageStop = new Image(getClass().getResourceAsStream("/stop.png"));
        ImageView viewStop = new ImageView(imageStop); 
        viewStop.setFitWidth(20);
        viewStop.setFitHeight(20);
        buttonStop = new Button("Stop", viewStop);
        buttonStop.setPrefSize(100, 20);
        buttonStop.setDisable(true);
        buttonStop.setOnAction((event) -> {
			String confirm = FxDialogs.showConfirm( "STOP button just being hit...", "Do you really want to stop the process?", "No", "Yes");
			if(confirm.equals("No")){
				return;
			}
			myass.setStopSignal(true);
        	buttonStop.setDisable(true);


        	//buttonRestart.setDisable(false);

        });
        
        
        hbox.getChildren().addAll(buttonStart, buttonStop);
        
        return hbox;
    }
    private boolean checkFileFromTextField(TextField tf) {
		String _path = tf.getText().trim();				
		if (_path.equals("")){
			FxDialogs.showWarning("File not found!", "Text field must not be empty!");
			tf.requestFocus();
			return false;
		}
		File _file = new File(_path);
		if (!_file.isFile()){
			FxDialogs.showWarning("File not found!", "File \"" + _path + "\" does not exist!");
			tf.requestFocus();
			return false;
		}
		return true;
    }
    private boolean checkFolderFromTextField(TextField tf) {
		String _path = tf.getText().trim();				
		if (_path.equals("")){
			FxDialogs.showWarning("Directory not found!", "Text field must not be empty!");
			tf.requestFocus();
			return false;
		}
		File _file = new File(_path);
		if (!_file.isDirectory()){
			FxDialogs.showWarning("File not found!", "Directory \"" + _path + "\" does not exist!");
			tf.requestFocus();
			return false;
		}
		return true;
    }
    private boolean checkingAndSetting() {	
    	if(myass==null)
    		myass=new HybridAssembler();
		//1. validate short-read assembly input
		if(!checkFileFromTextField(shortInputTF))
			return false;
		myass.setShortReadsInput(shortInputTF.getText());
		if(!checkFileFromTextField(longInputTF))
			return false;
		myass.setLongReadsInput(longInputTF.getText());
		if(!checkFolderFromTextField(outputTF))
			return false;
		myass.setPrefix(outputTF.getText());
		if(shortInputFormatCombo.getValue().equals("fasta/fastq")) {
			if(!checkFolderFromTextField(mm2PathTF))
				return false;
			if(!myass.checkMinimap2())
				return false;
			myass.setMinimapPath(mm2PathTF.getText());
			myass.setMinimapOpts(mm2OptTF.getText());
			
			Alignment.MIN_QUAL=Integer.valueOf(minQualTF.getText());
		}
		try{
			System.setProperty("usr.dir", myass.getPrefix());
		}
		catch(NullPointerException | IllegalArgumentException | SecurityException e ){
			e.printStackTrace();
			FxDialogs.showWarning("Illegal output folder!", "Please specify another output destination");
			outputTF.requestFocus();
			return false;
		}
		
		if(!myass.prepareShortReadsProcess(true)) {
			FxDialogs.showWarning("Warning", "Problems preparing assembly graph file. Check stderr!");
			return false;
		}

		if(!myass.prepareLongReadsProcess()) {
			FxDialogs.showWarning("Warning", "Problems preparing long-reads data. Check stderr");
			return false;
		}

    	GraphUtil.redrawGraphComponents(myass.simGraph);
    	myass.setReady(true);
    	return true;
    }
    
    
    private final int LeftPaneWidth=360;

    /*
     * Creates a VBox with a list of parameter settings
     */
    private VBox addVBox(Stage pStage, Stage gStage) {
        
        VBox vbox = new VBox();
        vbox.setPadding(new Insets(10)); // Set all sides to 10
        vbox.setSpacing(8);              // Gap between nodes
 
        final Text title = new Text("Settings");
        title.setFont(Font.font("Arial", FontWeight.BOLD, 15));
        vbox.getChildren().add(title);
        final Separator sep1 = new Separator();
        sep1.setMaxWidth(LeftPaneWidth);
        vbox.getChildren().add(1, sep1);
        
        vbox.getChildren().add(addInputPane(pStage, gStage));
        vbox.setSpacing(5);
        final Separator sep2 = new Separator();
        sep2.setMaxWidth(LeftPaneWidth);
        vbox.getChildren().add(3, sep2);

        
        vbox.getChildren().add(addOutputPane(pStage));
        vbox.setSpacing(5);
        final Separator sep3 = new Separator();
        sep3.setMaxWidth(LeftPaneWidth);
        vbox.getChildren().add(5, sep3);
        
        vbox.getChildren().add(addOptionPane(pStage));
               
        
        return vbox;
    }
    //TODO: put load graph button here to link with gStage.show()
    private GridPane addInputPane(Stage pStage, Stage gStage) {
    	GridPane inputPane = createFixGridPane(LeftPaneWidth, 5);
    	
    	final Label inputLabel = new Label("Input:");
    	inputLabel.setFont(Font.font("Roman", FontWeight.BOLD, 12));
    	inputLabel.setStyle("-fx-underline:true");
    	GridPane.setConstraints(inputLabel, 0,0);
    	inputPane.getChildren().add(inputLabel);
    	   		
    	final Label shortInputLabel = new Label("1. Pre-assemblies:");
    	shortInputLabel.setFont(Font.font("Roman", FontWeight.SEMI_BOLD, 12));
    	GridPane.setConstraints(shortInputLabel, 0,1,2,1);
    	inputPane.getChildren().add(shortInputLabel);
    	
    	shortInputFormatCombo=new ComboBox<String>();
        shortInputFormatCombo.getItems().addAll("fastg", "gfa");   
        shortInputFormatCombo.setValue(myass.getShortReadsInputFormat());
        shortInputFormatCombo.valueProperty().addListener((obs_val, old_val, new_val) -> {
        	myass.setShortReadsInputFormat(new_val);
        	shortInputTF.setText("");
        	myass.setShortReadsInput("");
        });
        GridPane.setConstraints(shortInputFormatCombo, 2, 1, 2, 1);
        inputPane.getChildren().add(shortInputFormatCombo);
    	
    	shortInputTF = new TextField("");
    	shortInputTF.setPromptText("Enter file name for assembly graph...");
    	if(!myass.getShortReadsInput().isEmpty())
    		shortInputTF.setText(myass.getShortReadsInput());
    	shortInputTF.setOnKeyPressed(e -> {
            if (e.getCode() == KeyCode.ENTER)  {
                buttonStart.requestFocus();
            }
    	});
    	//textField.setPrefWidth(250);
    	GridPane.setConstraints(shortInputTF, 0,2,4,1);
    	inputPane.getChildren().add(shortInputTF);
    	

    	shortInputBrowseButton = new ImageButton("/folder.png");
    	shortInputBrowseButton.setPrefSize(10, 10);
    	shortInputBrowseButton.setOnAction((event) -> {
       		FileChooser chooser = new FileChooser();
    		chooser.setTitle("Select assembly graph file");
    		File defaultFile = new File(shortInputTF.getText());
    		if(defaultFile.isFile())
    			chooser.setInitialFileName(defaultFile.getName());
    		if(defaultFile.getParentFile() !=null && defaultFile.getParentFile().isDirectory())
    			chooser.setInitialDirectory(defaultFile.getParentFile());
    		chooser.setSelectedExtensionFilter(
    				new ExtensionFilter("Assembly graph", 	"*.fastg", "*.FASTG", "*.gfa", "*.GFA"));
    		File selectedFile = chooser.showOpenDialog(pStage);
    		if(selectedFile != null){				
				try {
					myass.setShortReadsInput(selectedFile.getCanonicalPath());
					shortInputFormatCombo.setValue(myass.getShortReadsInputFormat());
					
					shortInputTF.setText(myass.getShortReadsInput());
				} catch (IOException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}	

    		}
        });
    	GridPane.setConstraints(shortInputBrowseButton, 4,2);
    	GridPane.setHalignment(shortInputBrowseButton,HPos.LEFT);
    	inputPane.getChildren().add(shortInputBrowseButton);
    	//inputPane.setGridLinesVisible(true);

    	graphCB = new CheckBox("Show assembly graph");
    	graphCB.setSelected(myass.simGraph!=null);
    	GridPane.setConstraints(graphCB, 0,3,5,1);
    	inputPane.getChildren().add(graphCB);
    	
    	final Label longInputLabel = new Label("2. Long-reads data:");
    	longInputLabel.setFont(Font.font("Roman", FontWeight.SEMI_BOLD, 12));
    	GridPane.setConstraints(longInputLabel, 0,4,2,1);
    	inputPane.getChildren().add(longInputLabel);
    	
    	longInputFormatCombo = new ComboBox<>();
        longInputFormatCombo.getItems().addAll("fasta/fastq", "sam/bam");   
        longInputFormatCombo.setValue(myass.getLongReadsInputFormat());
        longInputFormatCombo.valueProperty().addListener((obs_val, old_val, new_val) -> {
        	myass.setLongReadsInputFormat(new_val);
        	longInputTF.setText("");
        	myass.setLongReadsInput("");
        });
        GridPane.setConstraints(longInputFormatCombo, 2, 4, 2, 1);
        inputPane.getChildren().add(longInputFormatCombo);
    	
    	longInputTF = new TextField("");
    	longInputTF.setPromptText("Enter file name of long-reads data...");
    	if(!myass.getLongReadsInput().isEmpty())
    		longInputTF.setText(myass.getLongReadsInput());
    	longInputTF.setOnKeyPressed(e -> {
            if (e.getCode() == KeyCode.ENTER)  {
                buttonStart.requestFocus();
            }
    	});
    	GridPane.setConstraints(longInputTF, 0,5,4,1);
    	inputPane.getChildren().add(longInputTF);
    	
    	longReadsBrowseButton = new ImageButton("/folder.png");
    	longReadsBrowseButton.setPrefSize(10, 10);
    	longReadsBrowseButton.setDisable(!graphCB.isSelected());
    	longReadsBrowseButton.setOnAction((event) -> {
    		FileChooser chooser = new FileChooser();
    		chooser.setTitle("Select long-reads data file");
    		File defaultFile = new File(longInputTF.getText());
    		if(defaultFile.isFile())
    			chooser.setInitialFileName(defaultFile.getName());
    		if(defaultFile.getParentFile()!=null && defaultFile.getParentFile().isDirectory())
    			chooser.setInitialDirectory(defaultFile.getParentFile());
    		chooser.setSelectedExtensionFilter(
    				new ExtensionFilter("Long-reads data", 	"*.fastq", "*.fasta", "*.fq", "*.fa", "*.fna", "*.sam", "*.bam" , 
    														"*.FASTQ", "*.FASTA", "*.FQ", "*.FA", "*.FNA", "*.SAM", "*.BAM"));
    		File selectedFile = chooser.showOpenDialog(pStage);
    		if(selectedFile != null){
				try {
					myass.setLongReadsInput(selectedFile.getCanonicalPath());
					longInputFormatCombo.setValue(myass.getLongReadsInputFormat());
					longInputTF.setText(myass.getLongReadsInput());	
				} catch (IOException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}

    		}
        });
    	
    	
    	GridPane.setConstraints(longReadsBrowseButton, 4,5);
    	GridPane.setHalignment(longReadsBrowseButton,HPos.LEFT);
    	inputPane.getChildren().add(longReadsBrowseButton);
    	
    	graphCB.selectedProperty().addListener(
                (obs_val,old_val,new_val) -> {
                	tabPane.getSelectionModel().select(0);
                });	
    	
		return inputPane;
	}
    private GridPane addOutputPane(Stage stage) {
    	GridPane outputPane = createFixGridPane(LeftPaneWidth, 5);
    	
    	final Label outputLabel = new Label("Output:");
    	outputLabel.setFont(Font.font("Roman", FontWeight.BOLD, 12));
    	outputLabel.setStyle("-fx-underline:true");
    	GridPane.setConstraints(outputLabel, 0,0);
    	outputPane.getChildren().add(outputLabel);
    	
    	outputTF = new TextField("");
    	if(!myass.getPrefix().isEmpty())
    		outputTF.setText(myass.getPrefix());
    	outputTF.setPromptText("Enter name for output file...");
    	outputTF.setOnKeyPressed(e -> {
            if (e.getCode() == KeyCode.ENTER)  {
                buttonStart.requestFocus();
            }
    	});
    	GridPane.setConstraints(outputTF, 0,1,4,1);
    	outputPane.getChildren().add(outputTF);
    	

    	outputBrowseButton = new ImageButton("/folder.png");
    	outputBrowseButton.setPrefSize(10, 10);
//    	outputBrowseButton.setDisable(assembler.output.equals("-"));
    	outputBrowseButton.setOnAction((event) -> {
    		DirectoryChooser chooser = new DirectoryChooser();
    		chooser.setTitle("Save output to a destination folder");
    		File defaultDirectory=new File(outputTF.getText());
    		if(defaultDirectory.isDirectory())
    			chooser.setInitialDirectory(defaultDirectory);
    		File selectedDirectory=chooser.showDialog(stage);
    		if(selectedDirectory != null) {
    			myass.setPrefix(selectedDirectory.getPath());
    			outputTF.setText(myass.getPrefix());
    		}

        });
    	GridPane.setConstraints(outputBrowseButton, 4,1);
    	GridPane.setHalignment(outputBrowseButton, HPos.LEFT);
    	outputPane.getChildren().add(outputBrowseButton);
    	
    	overwriteCB = new CheckBox("Overwrite existing files if needed");
    	overwriteCB.selectedProperty().addListener(
    			(obs_val, old_val, new_val) -> {
    				myass.setOverwrite(new_val);
    			}		
		);
    	GridPane.setConstraints(overwriteCB, 0, 2, 4, 1);
    	outputPane.getChildren().add(overwriteCB);
    	
		return outputPane;
	}
    private GridPane addOptionPane(Stage stage) {
    	GridPane optionPane = createFixGridPane(LeftPaneWidth, 5);
    	
    	final Label optLabel = new Label("Alignment options:");
    	optLabel.setFont(Font.font("Roman", FontWeight.BOLD, 12));
    	optLabel.setStyle("-fx-underline:true");
    	
    	GridPane.setConstraints(optLabel, 0,0,4,1);
    	optionPane.getChildren().add(optLabel);
    	
    	final Label label1 = new Label("Path to ./minimap2: "),
					label2= new Label("Parameters setting:");
    	
    	GridPane.setConstraints(label1, 0,1,4,1);
    	optionPane.getChildren().add(label1);
    	
       	mm2PathTF = new TextField("");
       	mm2PathTF.setPromptText("Enter path to minimap2...");
       	if(!myass.getMinimapPath().isEmpty())
       		mm2PathTF.setText(myass.getMinimapPath());
       	mm2PathTF.setOnKeyPressed(e -> {
            if (e.getCode() == KeyCode.ENTER)  {
            	myass.setMinimapPath(mm2PathTF.getText());
                buttonStart.requestFocus();
            }
    	});
    	GridPane.setConstraints(mm2PathTF, 0,2,4,1);
    	optionPane.getChildren().add(mm2PathTF);
    	

    	mm2BrowseButton = new ImageButton("/folder.png");
    	mm2BrowseButton.setPrefSize(10, 10);
    	mm2BrowseButton.setOnAction((event) -> {
    		DirectoryChooser chooser = new DirectoryChooser();
    		chooser.setTitle("Folder containing minimap2");
    		File defaultDirectory=new File(mm2PathTF.getText());
    		if(defaultDirectory.isDirectory())
    			chooser.setInitialDirectory(defaultDirectory);
    		File selectedDirectory=chooser.showDialog(stage);
    		if(selectedDirectory != null) {
    			myass.setMinimapPath(selectedDirectory.getPath());
    			mm2PathTF.setText(myass.getMinimapPath());
    		}

        });
    	GridPane.setConstraints(mm2BrowseButton,4,2);
    	optionPane.getChildren().add(mm2BrowseButton);
    	
    	GridPane.setConstraints(label2, 0,3,4,1);
    	optionPane.getChildren().add(label2);
    	
       	mm2OptTF = new TextField("");
       	mm2OptTF.setPromptText("Enter options to minimap2...");
       	if(!myass.getMinimapOpts().isEmpty())
       		mm2OptTF.setText(myass.getMinimapOpts());
       	mm2OptTF.setOnKeyPressed(e -> {
            if (e.getCode() == KeyCode.ENTER)  {
            	myass.setMinimapPath(mm2OptTF.getText());
                buttonStart.requestFocus();
            }
    	});
    	GridPane.setConstraints(mm2OptTF, 0,4,4,1);
    	optionPane.getChildren().add(mm2OptTF);
    	
    	
    	final Label labelQual = new Label("Must have quality greater than ");
    	GridPane.setConstraints(labelQual, 0,6,3,1);
    	optionPane.getChildren().add(labelQual);
    	
    	minQualTF = new TextField("");
    	minQualTF.setPromptText("min.");
    	minQualTF.setText(Integer.toString(Alignment.MIN_QUAL));
    	minQualTF.setOnKeyPressed(e -> {
            if (e.getCode() == KeyCode.ENTER)  {
                buttonStart.requestFocus();
            }
    	});
    	GridPane.setConstraints(minQualTF, 3,6);
    	optionPane.getChildren().add(minQualTF);
    	
    	optionPane.disableProperty().bind(longInputFormatCombo.valueProperty().isEqualTo("sam/bam"));  
		return optionPane;
	}
    
    private GridPane createFixGridPane(int width, int ncols){
        GridPane gridpane = new GridPane();
        for (int i = 0; i < ncols; i++) {
            ColumnConstraints column = new ColumnConstraints(1.0*width/ncols);
            gridpane.getColumnConstraints().add(column);
        }
        gridpane.setPadding(new Insets(10, 10, 10, 10));
        gridpane.setVgap(5);
        gridpane.setHgap(5);
        return gridpane;
    }
    private GridPane createAutoresizeGridPane(int ncols, int nrows){
        GridPane gridpane = new GridPane();
        for (int i = 0; i < ncols; i++) {
            ColumnConstraints column = new ColumnConstraints();
            column.setPercentWidth(100/ncols);
            gridpane.getColumnConstraints().add(column);
        }
        for (int i = 0; i < nrows; i++) {
            RowConstraints row = new RowConstraints();
            row.setPercentHeight(100/nrows);
            gridpane.getRowConstraints().add(row);
        }
        gridpane.setPadding(new Insets(5, 5, 5, 5));
        gridpane.setVgap(5);
        gridpane.setHgap(5);
        return gridpane;
    }
     
    /*
     * Creates a grid for the center region with 2 columns and 2 rows
     */
    private GridPane addAssemblyStatsPane() {
    	
 	   GridPane mainGrid = createAutoresizeGridPane(1,2);
       mainGrid.setStyle("-fx-background-color: #C0C0C0;");
       mainGrid.setPadding(new Insets(5, 100, 5, 100));
       mainGrid.setVgap(5);
       mainGrid.setHgap(5);
       
//     GridPane.setConstraints(countPane, 0, 1);
//     mainGrid.getChildren().add(countPane);
     
//     mainGrid.setGridLinesVisible(true);
   		
   		return mainGrid;
    }
    
    /*
     * Create a Grid Pane for barcode analysis
     */
    private GridPane addGraphResolverPane(){  		
    	GridPane mainGrid = createAutoresizeGridPane(1,1);
		mainGrid.setStyle("-fx-background-color: #C0C0C0;");
		
		
		ThreadProxyPipe pipe = new ThreadProxyPipe() ;
		pipe.init(myass.simGraph);
		Viewer graphViewer = new FxViewer(pipe);
		System.setProperty("org.graphstream.ui", "javafx");

		FxDefaultView view = new FxDefaultView(graphViewer, "npGraph", new FxGraphRenderer());
		graphViewer.addView(view);
		graphViewer.enableAutoLayout();
		graphViewer.setCloseFramePolicy(CloseFramePolicy.CLOSE_VIEWER);
		
		mainGrid.getChildren().add(view);
		
		return mainGrid;
    }
    
    /******************************************************************************************
     * ** Here are variables and controls for all plots ***************************************
     ******************************************************************************************/
	final TextField txtCompReads= new TextField("0"), 
					txtTempReads= new TextField("0"), 
					txt2DReads= new TextField("0");
	final TextField txtPFiles= new TextField("0"), 
					txtFFiles= new TextField("0"), 
					txtTFiles= new TextField("0");

	//private static boolean stillRun = true;

	public static void interupt(Exception e){
		//stillRun = false;
//		assembler.wait = false;
		FxDialogs.showError("Unexpected errors happened!", e.getMessage());
	}
	
	private void updateData(){
		 
        final ScheduledExecutorService scheduler 
            = Executors.newScheduledThreadPool(1);
 
//        scheduler.scheduleAtFixedRate(
//                new Runnable(){         
//            		int lastIndexLengths = 0;//, lastIndexLengths2D = 0, lastIndexLengthsComp = 0, lastIndexLengthsTemp = 0;
//            		int lastIndexQual2D = 0, lastIndexQualComp = 0, lastIndexQualTemp = 0;
//                    @Override
//                    public void run() {
//                    	if(!assembler.wait)
//                    		scheduler.shutdown();
//                    	
//            			Second period = new Second();
//            			allReadsCount.add(period, assembler.twoDCount,"2D");
//            			allReadsCount.add(period, assembler.tempCount,"template");
//            			allReadsCount.add(period, assembler.compCount,"complement");
//
//            
//            			Demultiplexer myDmplx = assembler.dmplx;
//            			if(myDmplx!=null){
//            				for(int i=0;i<myDmplx.readCount.length;i++){
//            					demultiplexedStackReadsCount.add(period, myDmplx.readCount[i], myDmplx.barCodes.get(i).getName());
//            					
//            					demultiplexedBarReadsCount.setValue(myDmplx.readCount[i], "Count", myDmplx.barCodes.get(i).getName());
//            				}
//            				
//            			}
//            			
//            			txtTFiles.setText(assembler.getTotalFilesNumber()+"");	                
//            			txtPFiles.setText(assembler.getOKFilesNumber()+"");
//            			txtFFiles.setText(assembler.getSkippedFilesNumber()+"");
//            
//            			txt2DReads.setText(assembler.twoDCount+"");
//            			txtTempReads.setText(assembler.tempCount+"");
//            			txtCompReads.setText(assembler.compCount+"");
//
//            
//            			int currentIndex = assembler.lengths.size();
//            
//            			if (currentIndex > lastIndexLengths){			
//            				int index = histoLengthDataSet.getSeriesIndex("Read Length");
//            				for (int i = lastIndexLengths; i < currentIndex;i++)
//            					histoLengthDataSet.addSeries(index, assembler.lengths.get(i));
//            
//            				lastIndexLengths = currentIndex;
//            
//            				histoLengthDataSet.notifyChanged();
//            			}
//            
//            			currentIndex = assembler.qual2D.size();
//            			if (currentIndex > lastIndexQual2D){
//            				int index = histoQualDataSet.getSeriesIndex("2D");
//            				for (int i = lastIndexQual2D; i < currentIndex;i++)
//            					histoQualDataSet.addSeries(index, assembler.qual2D.get(i));
//            
//            				lastIndexQual2D = currentIndex;
//            				histoQualDataSet.notifyChanged();
//            			}
//            			currentIndex = assembler.qualTemp.size();
//            			if (currentIndex > lastIndexQualTemp){
//            				int index = histoQualDataSet.getSeriesIndex("template");
//            				for (int i = lastIndexQualTemp; i < currentIndex;i++)
//            					histoQualDataSet.addSeries(index, assembler.qualTemp.get(i));
//            
//            				lastIndexQualTemp = currentIndex;
//            				histoQualDataSet.notifyChanged();
//            			}
//            			currentIndex = assembler.qualComp.size();
//            			if (currentIndex > lastIndexQualComp){
//            				int index = histoQualDataSet.getSeriesIndex("complement");
//            				for (int i = lastIndexQualComp; i < currentIndex;i++)
//            					histoQualDataSet.addSeries(index, assembler.qualComp.get(i));
//            
//            				lastIndexQualComp = currentIndex;
//            				histoQualDataSet.notifyChanged();
//            			}          
//                        
//                    }
//                }, 
//                1, 
//                1, 
//                TimeUnit.SECONDS);     
	}


	public static void main(String[] args) {
		HybridAssembler hbAss = new HybridAssembler();
		
		hbAss.setShortReadsInput("/home/sonhoanghguyen/Projects/scaffolding/data/spades_v3.10/EcK12S-careful/assembly_graph.gfa");
		hbAss.setLongReadsInput("/home/sonhoanghguyen/Projects/scaffolding/data/Eck12_ONT.fasta");
		hbAss.setMinimapPath("/home/sonhoanghguyen/.usr/local/bin/"); 
		
		NPGraphFX.setAssembler(hbAss);
		Application.launch(NPGraphFX.class,args);
	}
	
}

