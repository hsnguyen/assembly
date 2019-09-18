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
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutorService;

import org.graphstream.graph.Node;
import org.graphstream.stream.thread.ThreadProxyPipe;
import org.graphstream.ui.fx_viewer.FxDefaultView;
import org.graphstream.ui.fx_viewer.FxViewer;
import org.graphstream.ui.javafx.FxGraphRenderer;
import org.graphstream.ui.view.Viewer;
import org.graphstream.ui.view.Viewer.CloseFramePolicy;
import org.rtassembly.npgraph.Alignment;
import org.rtassembly.npgraph.BDGraph;
import org.rtassembly.npgraph.HybridAssembler;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import japsa.util.FxDialogs;
import japsa.util.ImageButton;
import javafx.application.Application;
import javafx.application.Platform;
import javafx.geometry.HPos;
import javafx.geometry.Insets;
import javafx.scene.Scene;
import javafx.scene.control.Button;
import javafx.scene.control.CheckBox;
import javafx.scene.control.ComboBox;
import javafx.scene.control.Label;
import javafx.scene.control.Separator;
import javafx.scene.control.TextField;
import javafx.scene.image.Image;
import javafx.scene.image.ImageView;
import javafx.scene.input.KeyCode;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.ColumnConstraints;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.RowConstraints;
import javafx.scene.layout.VBox;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.scene.text.Text;

import javafx.scene.chart.XYChart;
import javafx.animation.AnimationTimer;
import javafx.scene.chart.AreaChart;
import javafx.scene.chart.NumberAxis;
import javafx.geometry.Side;

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
    	    // All the parameters setting to the left 
    	    graphInputPane = addGraphInputPane(primaryStage, graphStage);
    	    longInputPane = addLongReadInputPane(primaryStage);
    	    outputPane = addOutputPane(primaryStage);
    	    alignmentPane = addAlignmentOptionPane(primaryStage);
    	    controlPane = addControlPane();
    	    
    	    pBorder.setLeft(addVBox());
            
            // Here the main content    
            pBorder.setCenter(addAssemblyStatsBox());
            
            Scene pscene = new Scene(pBorder);
//            pscene.getStylesheets().add("/main.css");
            primaryStage.setScene(pscene);
            primaryStage.setTitle("npGraph Dashboard");
            primaryStage.setOnCloseRequest(e -> {
            	myass.terminateAlignmentProcess();
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
    private Button 	buttonStart, buttonStop, buttonGraph;
    private ComboBox<String> longInputFormatCombo;
    private GridPane graphInputPane, longInputPane, outputPane, alignmentPane, controlPane;
    
    
    private final int LeftPaneWidth=360;

    /*
     * Creates a VBox with a list of parameter settings
     */
    private VBox addVBox() {
        
        VBox vbox = new VBox();
        vbox.setPadding(new Insets(10)); // Set all sides to 10
        vbox.setSpacing(8);              // Gap between nodes
 
        final Text title = new Text("Settings");
        title.setFont(Font.font("Arial", FontWeight.BOLD, 15));
        vbox.getChildren().add(title);
        final Separator sep1 = new Separator();
        sep1.setMaxWidth(LeftPaneWidth);
        vbox.getChildren().add(1, sep1);
        
        vbox.getChildren().add(graphInputPane);
        vbox.setSpacing(5);
        final Separator sep2 = new Separator();
        sep2.setMaxWidth(LeftPaneWidth);
        vbox.getChildren().add(3, sep2);
        
        vbox.getChildren().add(longInputPane);
        vbox.setSpacing(5);
        final Separator sep3 = new Separator();
        sep2.setMaxWidth(LeftPaneWidth);
        vbox.getChildren().add(5, sep3);

        
        vbox.getChildren().add(outputPane);
        vbox.setSpacing(5);
        final Separator sep4 = new Separator();
        sep4.setMaxWidth(LeftPaneWidth);
        vbox.getChildren().add(7, sep4);
        
        vbox.getChildren().add(alignmentPane);
        vbox.setSpacing(5);
        final Separator sep5 = new Separator();
        sep4.setMaxWidth(LeftPaneWidth);
        vbox.getChildren().add(9, sep5);
               
        vbox.getChildren().add(controlPane);
        return vbox;
    }
    //TODO: put load graph button here to link with gStage.show()
    private GridPane addGraphInputPane(Stage pStage, Stage gStage) {
    	GridPane inputPane = createFixGridPane(LeftPaneWidth, 5);
    	
    	final Label inputLabel = new Label("Assembly graph:");
    	inputLabel.setFont(Font.font("Roman", FontWeight.BOLD, 12));
    	inputLabel.setStyle("-fx-underline:true");
    	GridPane.setConstraints(inputLabel, 0,0,2,1);
    	inputPane.getChildren().add(inputLabel);
    	
        TextField shortInputTF = new TextField("");
    	shortInputTF.setPromptText("Enter file name for assembly graph...");
    	shortInputTF.textProperty().bindBidirectional(myass.shortReadsInputProperty());
    	GridPane.setConstraints(shortInputTF, 0,1,4,1);
    	inputPane.getChildren().add(shortInputTF);
    	
    	ComboBox<String> shortInputFormatCombo=new ComboBox<String>();
        shortInputFormatCombo.getItems().addAll("fastg", "gfa");   
        shortInputFormatCombo.valueProperty().bindBidirectional(myass.shortReadsInputFormatProperty());

        GridPane.setConstraints(shortInputFormatCombo, 2, 0, 2, 1);
        inputPane.getChildren().add(shortInputFormatCombo);

    	Button shortInputBrowseButton = new ImageButton("/folder.png");
    	shortInputBrowseButton.setPrefSize(10, 10);
    	shortInputBrowseButton.setOnAction((event) -> {
       		FileChooser chooser = new FileChooser();
    		chooser.setTitle("Select assembly graph file");
    		File defaultFile = new File(myass.getShortReadsInput());
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
				} catch (IOException e1) {
	    			FxDialogs.showWarning("Warning", "Error loading graph file. Please try again!");
					e1.printStackTrace();
				}	

    		}
        });
    	GridPane.setConstraints(shortInputBrowseButton, 4,1);
    	GridPane.setHalignment(shortInputBrowseButton,HPos.LEFT);
    	inputPane.getChildren().add(shortInputBrowseButton);
    	//inputPane.setGridLinesVisible(true);

    	CheckBox binCB = new CheckBox("Use external binning information");
    	binCB.setSelected(false);
    	GridPane.setConstraints(binCB, 0,2,4,1);
    	inputPane.getChildren().add(binCB);
    	
    	TextField binInputTF = new TextField("");
    	binInputTF.setPromptText("Enter file name binning information...");
    	binInputTF.textProperty().bindBidirectional(myass.binReadsInputProperty());
    	binInputTF.disableProperty().bind(binCB.selectedProperty().not());
    	GridPane.setConstraints(binInputTF, 0,3,4,1);
    	inputPane.getChildren().add(binInputTF);
    	

    	Button binInputBrowseButton = new ImageButton("/folder.png");
    	binInputBrowseButton.setPrefSize(10, 10);
    	binInputBrowseButton.setOnAction((event) -> {
       		FileChooser chooser = new FileChooser();
    		chooser.setTitle("Select metaBAT bin file");
    		File defaultFile = new File(binInputTF.getText());
    		if(defaultFile.isFile())
    			chooser.setInitialFileName(defaultFile.getName());
    		if(defaultFile.getParentFile() !=null && defaultFile.getParentFile().isDirectory())
    			chooser.setInitialDirectory(defaultFile.getParentFile());
    		File selectedFile = chooser.showOpenDialog(pStage);
    		if(selectedFile != null){				
				try {
					myass.setBinReadsInput(selectedFile.getCanonicalPath());
				} catch (IOException e1) {
	    			FxDialogs.showWarning("Warning", "Error loading bin file. Please try again!");
	    			e1.printStackTrace();
				}	

    		}
        });
    	binInputBrowseButton.disableProperty().bind(binCB.selectedProperty().not());
    	GridPane.setConstraints(binInputBrowseButton, 4,3);
    	GridPane.setHalignment(binInputBrowseButton,HPos.LEFT);
    	inputPane.getChildren().add(binInputBrowseButton);   	   

    	CheckBox spadesCB = new CheckBox("Use SPAdes induced paths");
    	spadesCB.setSelected(myass.getUseSPAdesPath());
    	spadesCB.selectedProperty().bindBidirectional(myass.useSPAdesPathProperty());
    	GridPane.setConstraints(spadesCB, 0,4,4,1);
    	inputPane.getChildren().add(spadesCB);
    	
    	
    	CheckBox graphCB = new CheckBox("Show graph");
    	graphCB.setSelected(myass.simGraph!=null);
    	GridPane.setConstraints(graphCB, 3,5,2,1);
    	inputPane.getChildren().add(graphCB);
    	
        Image imageLoad = new Image(getClass().getResourceAsStream("/load.png"));
        ImageView viewLoad = new ImageView(imageLoad); 
        viewLoad.setFitWidth(20);
        viewLoad.setFitHeight(20);
        buttonGraph = new Button("Load", viewLoad);
        buttonGraph.setPrefSize(100, 20);
        buttonGraph.setOnAction((event) -> {
    		if(!myass.prepareShortReadsProcess()) { //true if using SPAdes path (not necessary)
    			FxDialogs.showWarning("Warning", myass.getCheckLog());
    			return;
    		}    		
        	if(graphCB.isSelected())
        		gStage.show();        	
			//TODO: unhide further controls
        	longInputPane.setDisable(false);
        	outputPane.setDisable(false);
        	controlPane.setDisable(false);
        	graphInputPane.setDisable(true);
        	updateData();
        	
		});
    	GridPane.setConstraints(buttonGraph, 0,5,3,1);
    	inputPane.getChildren().add(buttonGraph);
    	
		return inputPane;
	}
    
    private GridPane addLongReadInputPane(Stage pStage) {
    	GridPane inputPane = createFixGridPane(LeftPaneWidth, 5);
    	
    	final Label inputLabel = new Label("Long-read data:");
    	inputLabel.setFont(Font.font("Roman", FontWeight.BOLD, 12));
    	inputLabel.setStyle("-fx-underline:true");
    	GridPane.setConstraints(inputLabel, 0,0,2,1);
    	inputPane.getChildren().add(inputLabel);  	
    	
    	TextField longInputTF = new TextField("");
    	longInputTF.setPromptText("Enter file name of long-reads data...");
    	longInputTF.textProperty().bindBidirectional(myass.longReadsInputProperty());
    	GridPane.setConstraints(longInputTF, 0,1,4,1);
    	inputPane.getChildren().add(longInputTF);
    	
    	longInputFormatCombo = new ComboBox<>();
        longInputFormatCombo.getItems().addAll("fasta/fastq", "sam/bam");   
        longInputFormatCombo.valueProperty().bindBidirectional(myass.longReadsInputFormatProperty());
        GridPane.setConstraints(longInputFormatCombo, 2, 0, 2, 1);
        inputPane.getChildren().add(longInputFormatCombo);
    	
    	Button longReadsBrowseButton = new ImageButton("/folder.png");
    	longReadsBrowseButton.setPrefSize(10, 10);
    	longReadsBrowseButton.setOnAction((event) -> {
    		FileChooser chooser = new FileChooser();
    		chooser.setTitle("Select long-reads data file");
    		File defaultFile = new File(myass.getLongReadsInput());
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
				} catch (IOException e1) {
        			FxDialogs.showWarning("File not found!", "Please specify another long-read data file!");
					e1.printStackTrace();
				}

    		}
        });
    	
    	
    	GridPane.setConstraints(longReadsBrowseButton, 4,1);
    	GridPane.setHalignment(longReadsBrowseButton,HPos.LEFT);
    	inputPane.getChildren().add(longReadsBrowseButton);
    	
    	inputPane.setDisable(true);
		return inputPane;
	}
    
    private GridPane addOutputPane(Stage stage) {
    	GridPane outputPane = createFixGridPane(LeftPaneWidth, 5);
    	
    	final Label outputLabel = new Label("Output:");
    	outputLabel.setFont(Font.font("Roman", FontWeight.BOLD, 12));
    	outputLabel.setStyle("-fx-underline:true");
    	GridPane.setConstraints(outputLabel, 0,0);
    	outputPane.getChildren().add(outputLabel);
    	
    	TextField outputTF = new TextField("");
    	outputTF.setPromptText("Enter name for output file...");
    	outputTF.textProperty().bindBidirectional(myass.prefixProperty());
    	GridPane.setConstraints(outputTF, 0,1,4,1);
    	outputPane.getChildren().add(outputTF);
    	

    	Button outputBrowseButton = new ImageButton("/folder.png");
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
        		try{
        			System.setProperty("usr.dir", myass.getPrefix());
        		}
        		catch(NullPointerException | IllegalArgumentException | SecurityException exception ){
        			exception.printStackTrace();
        			FxDialogs.showWarning("Illegal output folder!", "Please specify another output destination");
        			return;
        		}
    		}

        });
    	GridPane.setConstraints(outputBrowseButton, 4,1);
    	GridPane.setHalignment(outputBrowseButton, HPos.LEFT);
    	outputPane.getChildren().add(outputBrowseButton);
    	
    	CheckBox overwriteCB = new CheckBox("Overwrite existing index files");
    	overwriteCB.selectedProperty().bindBidirectional(myass.overwriteProperty());
    	GridPane.setConstraints(overwriteCB, 0, 2, 4, 1);
    	outputPane.getChildren().add(overwriteCB);
    	
    	outputPane.setDisable(true);
		return outputPane;
	}
    private GridPane addAlignmentOptionPane(Stage stage) {
    	GridPane optionPane = createFixGridPane(LeftPaneWidth, 5);
    	
    	final Label optLabel = new Label("Alignment options:");
    	optLabel.setFont(Font.font("Roman", FontWeight.BOLD, 12));
    	optLabel.setStyle("-fx-underline:true");
    	
    	GridPane.setConstraints(optLabel, 0,0,4,1);
    	optionPane.getChildren().add(optLabel);
    	
    	        	
    	ComboBox<String> algCombo=new ComboBox<String>();
    	algCombo.getItems().addAll("minimap2", "bwa");   
    	algCombo.valueProperty().bindBidirectional(myass.alignerProperty());
        GridPane.setConstraints(algCombo, 2, 0, 2, 1);
        optionPane.getChildren().add(algCombo);         	
    	
    	
       	TextField algOptTF = new TextField("");
       	algOptTF.setPromptText("Enter options to aligner...");
       	algOptTF.textProperty().bindBidirectional(myass.alignerOptProperty());
    	GridPane.setConstraints(algOptTF, 0,1,4,1);
    	optionPane.getChildren().add(algOptTF);
    	
    	
    	final Label labelQual = new Label("Must have quality greater than ");
    	GridPane.setConstraints(labelQual, 0,2,3,1);
    	optionPane.getChildren().add(labelQual);
    	
    	TextField minQualTF = new TextField("");
//    	minQualTF.setPromptText("min.");
    	minQualTF.setText(Integer.toString(Alignment.MIN_QUAL));
    	minQualTF.setOnKeyPressed(e -> {
            if (e.getCode() == KeyCode.ENTER)  {
                buttonStart.requestFocus();
            }
    	});
    	GridPane.setConstraints(minQualTF, 3,2);
    	optionPane.getChildren().add(minQualTF);
    	
    	final Label consLabel = new Label("Consensus (optional):");
    	consLabel.setFont(Font.font("Roman", FontWeight.BOLD, 12));
    	consLabel.setStyle("-fx-underline:true");
    	GridPane.setConstraints(consLabel, 0,4,4,1);
    	optionPane.getChildren().add(consLabel);
    	
    	
    	ComboBox<String> consCombo=new ComboBox<String>();
    	consCombo.getItems().addAll("kalign", "none");   
    	consCombo.valueProperty().bindBidirectional(myass.msaProperty());
        GridPane.setConstraints(consCombo, 2, 4, 2, 1);
        optionPane.getChildren().add(consCombo); 
        
    	final Label labelMSA = new Label("Minium coverage for bridging ");
    	GridPane.setConstraints(labelMSA, 0,5,3,1);
    	optionPane.getChildren().add(labelMSA);
    	
    	TextField msaTF = new TextField("");
    	msaTF.setText(Integer.toString(BDGraph.MIN_SUPPORT));
    	msaTF.setOnKeyPressed(e -> {
            if (e.getCode() == KeyCode.ENTER)  {
                buttonStart.requestFocus();
            }
    	});
    	GridPane.setConstraints(msaTF, 3,5);
    	optionPane.getChildren().add(msaTF);
    	
    	
           	    	
    	optionPane.disableProperty().bind(	longInputPane.disabledProperty()
    										.or(longInputFormatCombo.valueProperty().isEqualTo("sam/bam")));  
		return optionPane;
	}
    
    private GridPane addControlPane() {
    	GridPane controlPane = createFixGridPane(LeftPaneWidth, 5);
    	
        Image imageStart = new Image(getClass().getResourceAsStream("/start.png"));
        ImageView viewStart = new ImageView(imageStart); 
        viewStart.setFitWidth(20);
        viewStart.setFitHeight(20);
        buttonStart = new Button("Start", viewStart);
        buttonStart.setPrefSize(100, 20);
        buttonStart.setOnAction((event) -> {
			String message=myass.getCheckLog();
    		if(!myass.prepareLongReadsProcess()){
    			FxDialogs.showError("Error", myass.getCheckLog());
    			return;
    		}else if(message!=null && !message.isEmpty())
    			FxDialogs.showWarning("Warning", myass.getCheckLog());

        	myass.setReady(true);
        	
			//Start running
			longInputPane.setDisable(true);
			outputPane.setDisable(true);
			
			buttonStart.setDisable(true);
			buttonStop.setDisable(false);
			
		});
        GridPane.setConstraints(buttonStart, 0,0,2,1);
    	controlPane.getChildren().add(buttonStart);
        
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
			if(!executor.isTerminated())
				executor.shutdown();
        	buttonStop.setDisable(true);

        	//buttonRestart.setDisable(false);

        });
        GridPane.setConstraints(buttonStop, 3,0,2,1);
    	controlPane.getChildren().add(buttonStop);
    	
    	controlPane.setDisable(true);
        return controlPane;
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
     
    
    /************************************************************************
     ********************** This is for the stats panel *********************
     ************************************************************************/
    private static final int MAX_DATA_POINTS = 500;
    
    private int xSeriesData = 0;
    private XYChart.Series<Number, Number> seriesN50 = new XYChart.Series<>();
    private XYChart.Series<Number, Number> seriesN75 = new XYChart.Series<>();
    private XYChart.Series<Number, Number> seriesMax = new XYChart.Series<>();
    private XYChart.Series<Number, Number> seriesNumCtgs = new XYChart.Series<>();
    private XYChart.Series<Number, Number> seriesNumCircularCtgs = new XYChart.Series<>();
    private ExecutorService executor;
    private ConcurrentLinkedQueue<Number> dataN50 = new ConcurrentLinkedQueue<>();
    private ConcurrentLinkedQueue<Number> dataN75 = new ConcurrentLinkedQueue<>();
    private ConcurrentLinkedQueue<Number> dataMax = new ConcurrentLinkedQueue<>();
    private ConcurrentLinkedQueue<Number> dataNumCtgs = new ConcurrentLinkedQueue<>();
    private ConcurrentLinkedQueue<Number> dataNumCircularCtgs = new ConcurrentLinkedQueue<>();

    private NumberAxis xAxis;
    
    @SuppressWarnings("unchecked")
	private VBox addAssemblyStatsBox(){
        VBox vbox = new VBox();
        vbox.setPadding(new Insets(20)); // Set all sides to 10
        vbox.setSpacing(8);              // Gap between nodes
        
        xAxis = new NumberAxis(0, MAX_DATA_POINTS, MAX_DATA_POINTS / 10);
        xAxis.setForceZeroInRange(false);
        xAxis.setAutoRanging(false);
        xAxis.setTickLabelsVisible(false);
        xAxis.setTickMarkVisible(false);
        xAxis.setMinorTickVisible(false);
        
        /*
         * Length stats chart
         */
        NumberAxis yAxis = new NumberAxis();
        yAxis.setTickLabelRotation(270);
        yAxis.setLabel("length (Kbp)");
        // Create a LineChart
        final AreaChart<Number, Number> lengthChart = new AreaChart<Number, Number>(xAxis, yAxis) {
            // Override to remove symbols on each data point
            @Override
            protected void dataItemAdded(Series<Number, Number> series, int itemIndex, Data<Number, Number> item) {
            }
        };
        
        lengthChart.setAnimated(false);
        lengthChart.setHorizontalGridLinesVisible(true);
        lengthChart.legendSideProperty().set(Side.TOP);

        // Set Name for Series
        seriesN50.setName("N50");
        seriesN75.setName("N75");
        seriesMax.setName("Max length");

        // Add Chart Series
        lengthChart.getData().addAll(seriesN50, seriesN75, seriesMax);
        
        vbox.getChildren().add(lengthChart);
		
		/*
		 * Number of contigs chart
		 */
		yAxis = new NumberAxis();
        yAxis.setTickLabelRotation(270);
        yAxis.setLabel("number of contigs");
         final AreaChart<Number, Number> numChart = new AreaChart<Number, Number>(xAxis, yAxis) {
             // Override to remove symbols on each data point
             @Override
             protected void dataItemAdded(Series<Number, Number> series, int itemIndex, Data<Number, Number> item) {
             }
         };       
         
         numChart.setAnimated(false);
         numChart.setHorizontalGridLinesVisible(true);
         
         seriesNumCtgs.setName("All contigs");
         seriesNumCircularCtgs.setName("Circular contigs");
         
         numChart.getData().addAll(seriesNumCtgs, seriesNumCircularCtgs);
         
         vbox.getChildren().add(numChart);
        
         lengthChart.getStylesheets().add("/chart1.css");
         numChart.getStylesheets().add("/chart2.css");

         return vbox;
    }

    /******************************************************************************************
     * ************ Create a Grid Pane for assembly graph *************************************
     ******************************************************************************************/
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
	//private static boolean stillRun = true;
    private class AddToQueue implements Runnable {
        public void run() {
            try {
                // add a item of data to queue
                dataN50.add(myass.observer.getN50()/1000); //because unit is Kbp
                dataN75.add(myass.observer.getN75()/1000);
                dataMax.add(myass.observer.getLongestContig()/1000);
                dataNumCtgs.add(myass.observer.getNumberOfSequences());
                dataNumCircularCtgs.add(myass.observer.getNumberOfCircularSequences());

                Thread.sleep(500);
                if(executor!=null && !executor.isShutdown())
                	executor.execute(this);
            } catch (InterruptedException ex) {
                ex.printStackTrace();
            }
        }
    }

    //-- Timeline gets called in the JavaFX Main thread
    private void prepareTimeline() {
        // Every frame to take any data from queue and add to chart
        new AnimationTimer() {
            @Override
            public void handle(long now) {
                addDataToSeries();
            }
        }.start();
    }

    private void addDataToSeries() {
        for (int i = 0; i < 20; i++) { //-- add 20 numbers to the plot+
            if (dataN50.isEmpty() || dataN75.isEmpty() || dataMax.isEmpty() 
        		|| dataNumCtgs.isEmpty() || dataNumCircularCtgs.isEmpty()) 
            	break;
            seriesN50.getData().add(new XYChart.Data<>(xSeriesData, dataN50.remove()));
            seriesN75.getData().add(new XYChart.Data<>(xSeriesData, dataN75.remove()));
            seriesMax.getData().add(new XYChart.Data<>(xSeriesData, dataMax.remove()));
            seriesNumCtgs.getData().add(new XYChart.Data<>(xSeriesData, dataNumCtgs.remove()));
            seriesNumCircularCtgs.getData().add(new XYChart.Data<>(xSeriesData++, dataNumCircularCtgs.remove()));
        }
        // remove points to keep us at no more than MAX_DATA_POINTS
        if (seriesN50.getData().size() > MAX_DATA_POINTS) {
        	seriesN50.getData().remove(0, seriesN50.getData().size() - MAX_DATA_POINTS);
        }
        if (seriesN75.getData().size() > MAX_DATA_POINTS) {
        	seriesN75.getData().remove(0, seriesN75.getData().size() - MAX_DATA_POINTS);
        }
        if (seriesMax.getData().size() > MAX_DATA_POINTS) {
        	seriesMax.getData().remove(0, seriesMax.getData().size() - MAX_DATA_POINTS);
        }
        if (seriesNumCtgs.getData().size() > MAX_DATA_POINTS) {
        	seriesNumCtgs.getData().remove(0, seriesNumCtgs.getData().size() - MAX_DATA_POINTS);
        }
        if (seriesNumCircularCtgs.getData().size() > MAX_DATA_POINTS) {
        	seriesNumCircularCtgs.getData().remove(0, seriesNumCircularCtgs.getData().size() - MAX_DATA_POINTS);
        }
        // update
        xAxis.setLowerBound(xSeriesData - MAX_DATA_POINTS);
        xAxis.setUpperBound(xSeriesData - 1);
    }
    
	
	private void updateData(){
        executor = Executors.newCachedThreadPool(new ThreadFactory() {
            @Override
            public Thread newThread(Runnable r) {
                Thread thread = new Thread(r);
                thread.setDaemon(true);
                return thread;
            }
        });

        AddToQueue addToQueue = new AddToQueue();
        executor.execute(addToQueue);
        //-- Prepare Timeline
        prepareTimeline();  
	}


	public static void main(String[] args) {
		HybridAssembler hbAss = new HybridAssembler();
		
////		//desktop IMB
////		hbAss.setShortReadsInput("/home/sonhoanghguyen/Projects/scaffolding/data/spades_v3.10/EcK12S-careful/assembly_graph.gfa");
////		hbAss.setLongReadsInput("/home/sonhoanghguyen/Projects/scaffolding/data/Eck12_ONT.fasta");
////		hbAss.setAlignerPath("/home/sonhoanghguyen/.usr/local/bin/"); 

		NPGraphFX.setAssembler(hbAss);
		Application.launch(NPGraphFX.class,args);
	}
	
}

