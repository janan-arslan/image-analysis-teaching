/**
 * Enhanced Cellpose Multi-Channel Segmentation with Colocalization Analysis
 * @author Janan Arslan
 * Part of the code includes Cellpose segmentation developed from Olivier Burri's template
 * 
 * Last edit: 11th June 2025
 * 
 * This script provides:
 * 1. Interactive channel selection for segmentation
 * 2. Multi-stain segmentation with combination analysis
 * 3. Comprehensive colocalization statistics
 * 4. Simplified, intuitive reporting
 */

import qupath.ext.biop.cellpose.Cellpose2D
import qupath.fx.dialogs.Dialogs
import qupath.lib.regions.RegionRequest
import ij.process.ByteProcessor
import ij.process.ImageProcessor
import java.awt.image.BufferedImage
import qupath.imagej.tools.IJTools
import qupath.lib.images.servers.ImageServer
import qupath.lib.objects.PathObject
import qupath.lib.images.PathImage
import qupath.lib.gui.dialogs.Dialogs as QuPathDialogs
import java.io.File
import java.io.FileWriter
import java.io.PrintWriter
import java.text.SimpleDateFormat
import java.util.Date

// Helper functions
def createObjectMask(PathImage pathImage, PathObject object, String objectType) {
    def bp = new ByteProcessor(pathImage.getImage().getWidth(), pathImage.getImage().getHeight())
    bp.setValue(1.0)

    if (objectType == "nucleus") {
        def roi = object.getNucleusROI()
        def roiIJ = IJTools.convertToIJRoi(roi, pathImage)
        bp.fill(roiIJ)
    } else if (objectType == "cytoplasm") {
        def nucleus = object.getNucleusROI()
        def roiIJNuc = IJTools.convertToIJRoi(nucleus, pathImage)
        def roi = object.getROI()
        def roiIJ = IJTools.convertToIJRoi(roi, pathImage)
        bp.fill(roiIJ)
        bp.setValue(0)
        bp.fill(roiIJNuc)
    } else {
        def roi = object.getROI()
        def roiIJ = IJTools.convertToIJRoi(roi, pathImage)
        bp.fill(roiIJ)
    }

    return bp
}

def calculateSummaryStats(values) {
    if (values.isEmpty()) return [mean: 0, median: 0, std: 0, min: 0, max: 0]
    
    values.sort()
    def mean = values.sum() / values.size()
    def median = values.size() % 2 == 0 ? 
                (values[values.size()/2 - 1] + values[values.size()/2]) / 2 : 
                values[values.size()/2]
    
    def variance = values.collect { (it - mean) ** 2 }.sum() / values.size()
    def std = Math.sqrt(variance)
    
    return [
        mean: mean,
        median: median,
        std: std,
        min: values[0],
        max: values[-1]
    ]
}

def saveHTMLReport(outputFolder, channelNames, segResults, comboResults, colocResults, pixelSize, diameter, model) {
    def reportFile = new File(outputFolder, "Colocalization_Analysis_Report.html")
    def writer = new PrintWriter(new FileWriter(reportFile))
    
    def timestamp = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss").format(new Date())
    def imageName = getCurrentImageData()?.getServer()?.getMetadata()?.getName() ?: "Unknown Image"
    
    writer.println("""<!DOCTYPE html>
<html>
<head>
    <title>Colocalization Analysis Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        .header { background-color: #f0f0f0; padding: 20px; border-radius: 5px; }
        .section { margin: 20px 0; }
        .stats-table { border-collapse: collapse; width: 100%; }
        .stats-table th, .stats-table td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        .stats-table th { background-color: #f2f2f2; }
        .result-box { background-color: #f9f9f9; padding: 15px; border-left: 4px solid #4CAF50; margin: 10px 0; }
        .warning { border-left-color: #ff9800; }
        .error { border-left-color: #f44336; }
    </style>
</head>
<body>
    <div class="header">
        <h1>Multi-Channel Colocalization Analysis Report</h1>
        <p><strong>Image:</strong> ${imageName}</p>
        <p><strong>Analysis Date:</strong> ${timestamp}</p>
        <p><strong>Cellpose Model:</strong> ${model} (diameter: ${diameter}px, pixel size: ${pixelSize}μm)</p>
    </div>
    
    <div class="section">
        <h2>Analysis Parameters</h2>
        <ul>
            <li>Channels analyzed: ${channelNames.join(', ')}</li>
            <li>Segmentation model: ${model}</li>
            <li>Expected diameter: ${diameter} pixels</li>
            <li>Pixel size: ${pixelSize} μm</li>
        </ul>
    </div>
    
    <div class="section">
        <h2>Segmentation Summary</h2>
        <table class="stats-table">
            <tr><th>Channel</th><th>Objects Detected</th></tr>""")
    
    channelNames.each { channel ->
        def count = segResults[channel]?.size() ?: 0
        writer.println("            <tr><td>${channel}</td><td>${count}</td></tr>")
    }
    
    writer.println("""        </table>
    </div>
    
    <div class="section">
        <h2>Combination Analysis</h2>
        <table class="stats-table">
            <tr><th>Combination</th><th>Co-positive Objects</th></tr>""")
    
    comboResults.each { combo, objects ->
        writer.println("            <tr><td>${combo}</td><td>${objects.size()}</td></tr>")
    }
    
    writer.println("""        </table>
    </div>
    
    <div class="section">
        <h2>Colocalization Statistics</h2>""")
    
    colocResults.each { analysisName, results ->
        def pearsonStats = calculateSummaryStats(results.pearson)
        def m1Stats = calculateSummaryStats(results.m1)
        def m2Stats = calculateSummaryStats(results.m2)
        
        def colocLevel = ""
        def boxClass = "result-box"
        if (pearsonStats.mean > 0.7) {
            colocLevel = "Strong positive"
            boxClass = "result-box"
        } else if (pearsonStats.mean > 0.3) {
            colocLevel = "Moderate positive"
            boxClass = "result-box"
        } else if (pearsonStats.mean > -0.3) {
            colocLevel = "Weak/No"
            boxClass = "result-box warning"
        } else {
            colocLevel = "Negative"
            boxClass = "result-box warning"
        }
        
        writer.println("""
        <div class="${boxClass}">
            <h3>${analysisName.replace('_', ' ').toUpperCase()} (n=${results.objectCount})</h3>
            <p><strong>Interpretation:</strong> ${colocLevel} colocalization</p>
            <table class="stats-table">
                <tr><th>Metric</th><th>Mean ± SD</th><th>Median</th><th>Range</th></tr>
                <tr><td>Pearson Correlation</td><td>${String.format('%.3f', (double)pearsonStats.mean)} ± ${String.format('%.3f', (double)pearsonStats.std)}</td><td>${String.format('%.3f', (double)pearsonStats.median)}</td><td>${String.format('%.3f', (double)pearsonStats.min)} - ${String.format('%.3f', (double)pearsonStats.max)}</td></tr>
                <tr><td>Manders M1</td><td>${String.format('%.3f', (double)m1Stats.mean)} ± ${String.format('%.3f', (double)m1Stats.std)}</td><td>${String.format('%.3f', (double)m1Stats.median)}</td><td>${String.format('%.3f', (double)m1Stats.min)} - ${String.format('%.3f', (double)m1Stats.max)}</td></tr>
                <tr><td>Manders M2</td><td>${String.format('%.3f', (double)m2Stats.mean)} ± ${String.format('%.3f', (double)m2Stats.std)}</td><td>${String.format('%.3f', (double)m2Stats.median)}</td><td>${String.format('%.3f', (double)m2Stats.min)} - ${String.format('%.3f', (double)m2Stats.max)}</td></tr>
            </table>
        </div>""")
    }
    
    writer.println("""
    </div>
    
    <div class="section">
        <h2>Analysis Notes</h2>
        <ul>
            <li>Pearson correlation coefficient ranges from -1 (perfect negative correlation) to +1 (perfect positive correlation)</li>
            <li>Manders M1: Fraction of channel 1 intensity that colocalizes with channel 2</li>
            <li>Manders M2: Fraction of channel 2 intensity that colocalizes with channel 1</li>
            <li>Strong colocalization: Pearson > 0.7</li>
            <li>Moderate colocalization: Pearson 0.3-0.7</li>
            <li>Weak/No colocalization: Pearson -0.3 to 0.3</li>
        </ul>
    </div>
    
</body>
</html>""")
    
    writer.close()
    println("✓ HTML report saved: ${reportFile.getName()}")
}

def saveCSVResults(outputFolder, colocResults) {
    colocResults.each { analysisName, results ->
        def csvFile = new File(outputFolder, "${analysisName}_detailed_results.csv")
        def writer = new PrintWriter(new FileWriter(csvFile))
        
        writer.println("Object_ID,Pearson_Correlation,Manders_M1,Manders_M2,Intensity_Ratio,Ch1_Mean_Intensity,Ch2_Mean_Intensity")
        
        for (int i = 0; i < results.pearson.size(); i++) {
            writer.println("${i+1},${results.pearson[i]},${results.m1[i]},${results.m2[i]},${results.intensityRatio[i]},${results.ch1Intensity[i]},${results.ch2Intensity[i]}")
        }
        
        writer.close()
        println("✓ CSV results saved: ${csvFile.getName()}")
    }
}

def saveObjectMeasurements(outputFolder, colocResults) {
    def csvFile = new File(outputFolder, "Object_Measurements_Summary.csv")
    def writer = new PrintWriter(new FileWriter(csvFile))
    
    writer.println("Analysis,Object_Count,Pearson_Mean,Pearson_Std,M1_Mean,M1_Std,M2_Mean,M2_Std")
    
    colocResults.each { analysisName, results ->
        def pearsonStats = calculateSummaryStats(results.pearson)
        def m1Stats = calculateSummaryStats(results.m1)
        def m2Stats = calculateSummaryStats(results.m2)
        
        writer.println("${analysisName},${results.objectCount},${String.format('%.4f', (double)pearsonStats.mean)},${String.format('%.4f', (double)pearsonStats.std)},${String.format('%.4f', (double)m1Stats.mean)},${String.format('%.4f', (double)m1Stats.std)},${String.format('%.4f', (double)m2Stats.mean)},${String.format('%.4f', (double)m2Stats.std)}")
    }
    
    writer.close()
    println("✓ Object measurements saved: ${csvFile.getName()}")
}

def saveSummaryStatistics(outputFolder, channelNames, segResults, comboResults, colocResults) {
    def summaryFile = new File(outputFolder, "Analysis_Summary.txt")
    def writer = new PrintWriter(new FileWriter(summaryFile))
    
    def timestamp = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss").format(new Date())
    def imageName = getCurrentImageData()?.getServer()?.getMetadata()?.getName() ?: "Unknown Image"
    
    writer.println("MULTI-CHANNEL COLOCALIZATION ANALYSIS SUMMARY")
    writer.println("=" * 50)
    writer.println("Image: ${imageName}")
    writer.println("Analysis Date: ${timestamp}")
    writer.println()
    
    writer.println("SEGMENTATION RESULTS:")
    writer.println("-" * 30)
    channelNames.each { channel ->
        def count = segResults[channel]?.size() ?: 0
        writer.println("${channel}: ${count} objects")
    }
    writer.println()
    
    writer.println("COMBINATION ANALYSIS:")
    writer.println("-" * 30)
    comboResults.each { combo, objects ->
        writer.println("${combo}: ${objects.size()} co-positive objects")
    }
    writer.println()
    
    writer.println("COLOCALIZATION SUMMARY:")
    writer.println("-" * 30)
    colocResults.each { analysisName, results ->
        def pearsonStats = calculateSummaryStats(results.pearson)
        def colocLevel = ""
        if (pearsonStats.mean > 0.7) colocLevel = "Strong positive"
        else if (pearsonStats.mean > 0.3) colocLevel = "Moderate positive"
        else if (pearsonStats.mean > -0.3) colocLevel = "Weak/No"
        else colocLevel = "Negative"
        
        writer.println("${analysisName}: ${colocLevel} colocalization (Pearson: ${String.format('%.3f', (double)pearsonStats.mean)} ± ${String.format('%.3f', (double)pearsonStats.std)})")
    }
    
    writer.close()
    println("✓ Summary statistics saved: ${summaryFile.getName()}")
}

def performSingleCellAnalysis(outputFolder, channelNames, segResults, colocResults, selectedChannelNumbers) {
    println("\n=== PERFORMING SINGLE-CELL COLOCALIZATION ANALYSIS ===")
    
    // Create single-cell subfolder
    def singleCellFolder = new File(outputFolder, "Single_Cell_Analysis")
    singleCellFolder.mkdirs()
    
    // For each analysis pair, perform detailed single-cell analysis
    colocResults.each { analysisName, results ->
        println("Analyzing single cells for: ${analysisName}")
        
        // Extract channel indices and names from analysis name
        def parts = analysisName.split("_vs_")
        if (parts.size() != 2) return
        
        def ch1Name = parts[0]
        def ch2Name = parts[1]
        
        // Find channel indices
        def ch1Index = -1
        def ch2Index = -1
        channelNames.eachWithIndex { name, index ->
            if (name == ch1Name) ch1Index = selectedChannelNumbers[index]
            if (name == ch2Name) ch2Index = selectedChannelNumbers[index]
        }
        
        if (ch1Index == -1 || ch2Index == -1) {
            println("  ⚠️  Could not find channel indices for ${analysisName}")
            return
        }
        
        // Get objects for this analysis
        def analysisObjects = segResults[ch1Name] + segResults[ch2Name]
        analysisObjects = analysisObjects.unique()
        
        if (analysisObjects.isEmpty()) {
            println("  ⚠️  No objects found for ${analysisName}")
            return
        }
        
        // Perform detailed single-cell analysis
        def singleCellData = []
        def server = getCurrentImageData().getServer()
        
        analysisObjects.eachWithIndex { object, objIndex ->
            try {
                def roi = object.getROI()
                def request = RegionRequest.createInstance(server.getPath(), 1.0, roi)
                def pathImage = IJTools.convertToImagePlus(server, request)
                def imp = pathImage.getImage()
                
                def ch1Image = imp.getStack().getProcessor(ch1Index).convertToFloatProcessor()
                def ch2Image = imp.getStack().getProcessor(ch2Index).convertToFloatProcessor()
                def ch1Pixels = ch1Image.getPixels()
                def ch2Pixels = ch2Image.getPixels()
                
                def mask = createObjectMask(pathImage, object, "cell").getPixels()
                
                def ch1Values = []
                def ch2Values = []
                for (int i = 0; i < ch1Pixels.size(); i++) {
                    if (mask[i]) {
                        ch1Values.add(ch1Pixels[i])
                        ch2Values.add(ch2Pixels[i])
                    }
                }
                
                if (ch1Values.size() == 0) return
                
                // Calculate comprehensive single-cell metrics
                def ch1Mean = ch1Values.sum() / ch1Values.size()
                def ch2Mean = ch2Values.sum() / ch2Values.size()
                def ch1Max = ch1Values.max()
                def ch2Max = ch2Values.max()
                def ch1Std = Math.sqrt(ch1Values.collect { (it - ch1Mean) ** 2 }.sum() / ch1Values.size())
                def ch2Std = Math.sqrt(ch2Values.collect { (it - ch2Mean) ** 2 }.sum() / ch2Values.size())
                
                // Pearson correlation
                def pearsonNum = 0
                def ch1SumSq = 0
                def ch2SumSq = 0
                
                ch1Values.eachWithIndex { val1, i ->
                    def val2 = ch2Values[i]
                    def diff1 = val1 - ch1Mean
                    def diff2 = val2 - ch2Mean
                    pearsonNum += diff1 * diff2
                    ch1SumSq += diff1 * diff1
                    ch2SumSq += diff2 * diff2
                }
                
                def pearson = 0
                if (ch1SumSq > 0 && ch2SumSq > 0) {
                    pearson = pearsonNum / Math.sqrt(ch1SumSq * ch2SumSq)
                }
                
                // Manders coefficients with multiple thresholds
                def ch1Background = ch1Mean * 0.1  // Adaptive threshold
                def ch2Background = ch2Mean * 0.1
                
                def m1Num = 0, m1Den = 0, m2Num = 0, m2Den = 0
                def overlapPixels = 0
                def ch1PositivePixels = 0
                def ch2PositivePixels = 0
                
                ch1Values.eachWithIndex { val1, i ->
                    def val2 = ch2Values[i]
                    if (val1 > ch1Background) ch1PositivePixels++
                    if (val2 > ch2Background) ch2PositivePixels++
                    if (val1 > ch1Background && val2 > ch2Background) overlapPixels++
                    
                    if (val2 > ch2Background) m1Num += Math.max(val1 - ch1Background, 0)
                    if (val1 > ch1Background) m2Num += Math.max(val2 - ch2Background, 0)
                    m1Den += Math.max(val1 - ch1Background, 0)
                    m2Den += Math.max(val2 - ch2Background, 0)
                }
                
                def m1 = m1Den > 0 ? m1Num / m1Den : 0
                def m2 = m2Den > 0 ? m2Num / m2Den : 0
                def overlapFraction = ch1Values.size() > 0 ? overlapPixels / ch1Values.size() : 0
                def ch1Fraction = ch1Values.size() > 0 ? ch1PositivePixels / ch1Values.size() : 0
                def ch2Fraction = ch1Values.size() > 0 ? ch2PositivePixels / ch1Values.size() : 0
                
                // Additional colocalization metrics
                def intensityCorrelation = ch1Mean > 0 ? ch2Mean / ch1Mean : 0
                def colocalizationStrength = pearson * overlapFraction  // Combined metric
                
                // Spatial metrics
                def cellArea = roi.getArea()
                def cellPerimeter = roi.getLength()
                def cellCircularity = cellPerimeter > 0 ? 4 * Math.PI * cellArea / (cellPerimeter * cellPerimeter) : 0
                
                // Store comprehensive single-cell data
                singleCellData.add([
                    cellId: objIndex + 1,
                    cellArea: cellArea,
                    cellPerimeter: cellPerimeter,
                    cellCircularity: cellCircularity,
                    totalPixels: ch1Values.size(),
                    ch1Mean: ch1Mean,
                    ch1Max: ch1Max,
                    ch1Std: ch1Std,
                    ch1PositivePixels: ch1PositivePixels,
                    ch1Fraction: ch1Fraction,
                    ch2Mean: ch2Mean,
                    ch2Max: ch2Max,
                    ch2Std: ch2Std,
                    ch2PositivePixels: ch2PositivePixels,
                    ch2Fraction: ch2Fraction,
                    pearsonCorrelation: pearson,
                    mandersM1: m1,
                    mandersM2: m2,
                    overlapPixels: overlapPixels,
                    overlapFraction: overlapFraction,
                    intensityCorrelation: intensityCorrelation,
                    colocalizationStrength: colocalizationStrength,
                    ch1Background: ch1Background,
                    ch2Background: ch2Background
                ])
                
            } catch (Exception e) {
                println("  ⚠️  Error analyzing cell ${objIndex + 1}: ${e.getMessage()}")
            }
        }
        
        // Save single-cell detailed CSV
        if (!singleCellData.isEmpty()) {
            saveSingleCellCSV(singleCellFolder, analysisName, singleCellData, ch1Name, ch2Name)
            saveSingleCellSummary(singleCellFolder, analysisName, singleCellData, ch1Name, ch2Name)
            generateSingleCellReport(singleCellFolder, analysisName, singleCellData, ch1Name, ch2Name)
        }
        
        println("  ✓ Single-cell analysis complete: ${singleCellData.size()} cells analyzed")
    }
    
    println("✓ Single-cell analysis saved to: ${singleCellFolder.getName()}")
}

def saveSingleCellCSV(folder, analysisName, cellData, ch1Name, ch2Name) {
    def csvFile = new File(folder, "${analysisName}_SingleCell_Detailed.csv")
    def writer = new PrintWriter(new FileWriter(csvFile))
    
    // Write header
    writer.println("Cell_ID,Cell_Area,Cell_Perimeter,Cell_Circularity,Total_Pixels," +
                   "${ch1Name}_Mean,${ch1Name}_Max,${ch1Name}_Std,${ch1Name}_PositivePixels,${ch1Name}_Fraction," +
                   "${ch2Name}_Mean,${ch2Name}_Max,${ch2Name}_Std,${ch2Name}_PositivePixels,${ch2Name}_Fraction," +
                   "Pearson_Correlation,Manders_M1,Manders_M2,Overlap_Pixels,Overlap_Fraction," +
                   "Intensity_Correlation,Colocalization_Strength,${ch1Name}_Background,${ch2Name}_Background")
    
    // Write data
    cellData.each { cell ->
        writer.println("${cell.cellId},${String.format('%.2f', (double)cell.cellArea)},${String.format('%.2f', (double)cell.cellPerimeter)}," +
                      "${String.format('%.3f', (double)cell.cellCircularity)},${cell.totalPixels}," +
                      "${String.format('%.2f', (double)cell.ch1Mean)},${String.format('%.2f', (double)cell.ch1Max)},${String.format('%.2f', (double)cell.ch1Std)}," +
                      "${cell.ch1PositivePixels},${String.format('%.3f', (double)cell.ch1Fraction)}," +
                      "${String.format('%.2f', (double)cell.ch2Mean)},${String.format('%.2f', (double)cell.ch2Max)},${String.format('%.2f', (double)cell.ch2Std)}," +
                      "${cell.ch2PositivePixels},${String.format('%.3f', (double)cell.ch2Fraction)}," +
                      "${String.format('%.4f', (double)cell.pearsonCorrelation)},${String.format('%.4f', (double)cell.mandersM1)},${String.format('%.4f', (double)cell.mandersM2)}," +
                      "${cell.overlapPixels},${String.format('%.3f', (double)cell.overlapFraction)}," +
                      "${String.format('%.3f', (double)cell.intensityCorrelation)},${String.format('%.4f', (double)cell.colocalizationStrength)}," +
                      "${String.format('%.2f', (double)cell.ch1Background)},${String.format('%.2f', (double)cell.ch2Background)}")
    }
    
    writer.close()
    println("  ✓ Single-cell CSV saved: ${csvFile.getName()}")
}

def saveSingleCellSummary(folder, analysisName, cellData, ch1Name, ch2Name) {
    def summaryFile = new File(folder, "${analysisName}_SingleCell_Summary.txt")
    def writer = new PrintWriter(new FileWriter(summaryFile))
    
    writer.println("SINGLE-CELL COLOCALIZATION ANALYSIS SUMMARY")
    writer.println("=" * 60)
    writer.println("Analysis: ${analysisName}")
    writer.println("Channels: ${ch1Name} vs ${ch2Name}")
    writer.println("Cells analyzed: ${cellData.size()}")
    writer.println("Date: ${new SimpleDateFormat('yyyy-MM-dd HH:mm:ss').format(new Date())}")
    writer.println()
    
    // Calculate summary statistics for each metric
    def metrics = [
        "Pearson Correlation": cellData.collect { it.pearsonCorrelation },
        "Manders M1": cellData.collect { it.mandersM1 },
        "Manders M2": cellData.collect { it.mandersM2 },
        "Overlap Fraction": cellData.collect { it.overlapFraction },
        "Colocalization Strength": cellData.collect { it.colocalizationStrength },
        "${ch1Name} Mean Intensity": cellData.collect { it.ch1Mean },
        "${ch2Name} Mean Intensity": cellData.collect { it.ch2Mean },
        "Cell Area": cellData.collect { it.cellArea }
    ]
    
    writer.println("SINGLE-CELL STATISTICS:")
    writer.println("-" * 40)
    metrics.each { metricName, values ->
        def stats = calculateSummaryStats(values)
        writer.println("${metricName}:")
        writer.println("  Mean ± SD: ${String.format('%.4f', (double)stats.mean)} ± ${String.format('%.4f', (double)stats.std)}")
        writer.println("  Median [Range]: ${String.format('%.4f', (double)stats.median)} [${String.format('%.4f', (double)stats.min)} - ${String.format('%.4f', (double)stats.max)}]")
        writer.println()
    }
    
    // Cell population analysis
    writer.println("CELL POPULATION ANALYSIS:")
    writer.println("-" * 40)
    
    def strongColoc = cellData.findAll { it.pearsonCorrelation > 0.7 }.size()
    def moderateColoc = cellData.findAll { it.pearsonCorrelation > 0.3 && it.pearsonCorrelation <= 0.7 }.size()
    def weakColoc = cellData.findAll { it.pearsonCorrelation >= -0.3 && it.pearsonCorrelation <= 0.3 }.size()
    def negativeColoc = cellData.findAll { it.pearsonCorrelation < -0.3 }.size()
    
    writer.println("Strong colocalization (r > 0.7): ${strongColoc} cells (${String.format('%.1f', strongColoc * 100.0 / cellData.size())}%)")
    writer.println("Moderate colocalization (0.3 < r ≤ 0.7): ${moderateColoc} cells (${String.format('%.1f', moderateColoc * 100.0 / cellData.size())}%)")
    writer.println("Weak/No colocalization (-0.3 ≤ r ≤ 0.3): ${weakColoc} cells (${String.format('%.1f', weakColoc * 100.0 / cellData.size())}%)")
    writer.println("Negative colocalization (r < -0.3): ${negativeColoc} cells (${String.format('%.1f', negativeColoc * 100.0 / cellData.size())}%)")
    
    writer.close()
    println("  ✓ Single-cell summary saved: ${summaryFile.getName()}")
}

def generateSingleCellReport(folder, analysisName, cellData, ch1Name, ch2Name) {
    def reportFile = new File(folder, "${analysisName}_SingleCell_Report.html")
    def writer = new PrintWriter(new FileWriter(reportFile))
    
    def timestamp = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss").format(new Date())
    
    writer.println("""<!DOCTYPE html>
<html>
<head>
    <title>Single-Cell Colocalization Report: ${analysisName}</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        .header { background-color: #e3f2fd; padding: 20px; border-radius: 5px; }
        .section { margin: 20px 0; }
        .stats-table { border-collapse: collapse; width: 100%; }
        .stats-table th, .stats-table td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        .stats-table th { background-color: #f2f2f2; }
        .highlight { background-color: #fff3e0; padding: 15px; border-left: 4px solid #ff9800; margin: 10px 0; }
        .cell-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; }
        .cell-card { border: 1px solid #ddd; padding: 15px; border-radius: 5px; }
        .high-coloc { border-left: 4px solid #4CAF50; }
        .med-coloc { border-left: 4px solid #FF9800; }
        .low-coloc { border-left: 4px solid #F44336; }
    </style>
</head>
<body>
    <div class="header">
        <h1>Single-Cell Colocalization Analysis</h1>
        <h2>${analysisName.replace('_', ' ')}</h2>
        <p><strong>Channels:</strong> ${ch1Name} vs ${ch2Name}</p>
        <p><strong>Cells Analyzed:</strong> ${cellData.size()}</p>
        <p><strong>Analysis Date:</strong> ${timestamp}</p>
    </div>
    
    <div class="section">
        <h2>Population Summary</h2>""")
    
    def strongColoc = cellData.findAll { it.pearsonCorrelation > 0.7 }.size()
    def moderateColoc = cellData.findAll { it.pearsonCorrelation > 0.3 && it.pearsonCorrelation <= 0.7 }.size()
    def weakColoc = cellData.findAll { it.pearsonCorrelation >= -0.3 && it.pearsonCorrelation <= 0.3 }.size()
    def negativeColoc = cellData.findAll { it.pearsonCorrelation < -0.3 }.size()
    
    writer.println("""
        <table class="stats-table">
            <tr><th>Colocalization Level</th><th>Cell Count</th><th>Percentage</th></tr>
            <tr><td>Strong (r > 0.7)</td><td>${strongColoc}</td><td>${String.format('%.1f', strongColoc * 100.0 / cellData.size())}%</td></tr>
            <tr><td>Moderate (0.3 < r ≤ 0.7)</td><td>${moderateColoc}</td><td>${String.format('%.1f', moderateColoc * 100.0 / cellData.size())}%</td></tr>
            <tr><td>Weak/None (-0.3 ≤ r ≤ 0.3)</td><td>${weakColoc}</td><td>${String.format('%.1f', weakColoc * 100.0 / cellData.size())}%</td></tr>
            <tr><td>Negative (r < -0.3)</td><td>${negativeColoc}</td><td>${String.format('%.1f', negativeColoc * 100.0 / cellData.size())}%</td></tr>
        </table>
    </div>
    
    <div class="section">
        <h2>Top Colocalizing Cells</h2>
        <div class="cell-grid">""")
    
    // Show top 10 colocalizing cells
    def topCells = cellData.sort { -it.pearsonCorrelation }.take(10)
    topCells.each { cell ->
        def colocClass = cell.pearsonCorrelation > 0.7 ? "high-coloc" : 
                        cell.pearsonCorrelation > 0.3 ? "med-coloc" : "low-coloc"
        writer.println("""
            <div class="cell-card ${colocClass}">
                <h4>Cell ${cell.cellId}</h4>
                <p><strong>Pearson r:</strong> ${String.format('%.3f', (double)cell.pearsonCorrelation)}</p>
                <p><strong>M1:</strong> ${String.format('%.3f', (double)cell.mandersM1)}, <strong>M2:</strong> ${String.format('%.3f', (double)cell.mandersM2)}</p>
                <p><strong>Overlap:</strong> ${String.format('%.1f', (double)cell.overlapFraction * 100)}%</p>
                <p><strong>Area:</strong> ${String.format('%.1f', (double)cell.cellArea)} px²</p>
            </div>""")
    }
    
    writer.println("""
        </div>
    </div>
    
    <div class="highlight">
        <h3>Key Insights</h3>
        <ul>
            <li>Cell population shows ${strongColoc > cellData.size() * 0.3 ? 'high' : moderateColoc > cellData.size() * 0.3 ? 'moderate' : 'low'} heterogeneity in colocalization</li>
            <li>Average colocalization strength: ${String.format('%.3f', (double)(cellData.collect{it.colocalizationStrength}.sum() / cellData.size()))}</li>
            <li>Most colocalizing cell: Cell ${topCells[0].cellId} (r = ${String.format('%.3f', (double)topCells[0].pearsonCorrelation)})</li>
        </ul>
    </div>
    
</body>
</html>""")
    
    writer.close()
    println("  ✓ Single-cell HTML report saved: ${reportFile.getName()}")
}

// Main analysis
def imageData = getCurrentImageData()
def server = imageData.getServer()

// Get available channels
def viewer = getCurrentViewer()
def display = viewer.getImageDisplay()
def availableChannels = display.availableChannels()

println("=== AVAILABLE CHANNELS ===")
def channelInfo = []
for (int i = 0; i < availableChannels.size(); i++) {
    def channelName = availableChannels[i].getName()
    channelInfo.add("${i+1}: ${channelName}")
    println("Channel ${i+1}: ${channelName}")
}

// User input for channel selection
def channelSelectionDialog = """
Available channels:
${channelInfo.join('\n')}

Please enter the channel numbers you want to analyze (comma-separated):
Example: 1,2,3 for channels 1, 2, and 3
"""

def selectedChannelsInput = Dialogs.showInputDialog("Channel Selection", channelSelectionDialog, "")
if (selectedChannelsInput == null || selectedChannelsInput.trim().isEmpty()) {
    println("No channels selected. Exiting.")
    return
}

// Parse selected channels
def selectedChannelNumbers = []
def selectedChannelNames = []
def selectedChannelNamesClean = []
try {
    selectedChannelsInput.split(',').each { channelStr ->
        def channelNum = Integer.parseInt(channelStr.trim())
        if (channelNum >= 1 && channelNum <= availableChannels.size()) {
            selectedChannelNumbers.add(channelNum)
            def fullChannelName = availableChannels[channelNum-1].getName()
            selectedChannelNames.add(fullChannelName)
            // Extract clean channel name (remove parentheses part)
            def cleanName = fullChannelName.replaceAll(/\s*\([^)]*\)$/, '').trim()
            selectedChannelNamesClean.add(cleanName)
        }
    }
} catch (Exception e) {
    Dialogs.showErrorMessage("Error", "Invalid channel input format. Please use comma-separated numbers.")
    return
}

if (selectedChannelNumbers.isEmpty()) {
    Dialogs.showErrorMessage("Error", "No valid channels selected.")
    return
}

println("\n=== SELECTED CHANNELS FOR ANALYSIS ===")
selectedChannelNames.eachWithIndex { name, index ->
    println("Channel ${selectedChannelNumbers[index]}: ${name} (Clean name: ${selectedChannelNamesClean[index]})")
}

// Output folder selection

def outputFolder = null
def useOutputFolder = Dialogs.showYesNoDialog("Output Folder", 
    "Would you like to specify an output folder to save analysis results?\n\n" +
    "This will save:\n" +
    "• Detailed analysis report (HTML)\n" +
    "• CSV data files\n" +
    "• Summary statistics\n" +
    "• Segmentation overlay images (optional)\n\n" +
    "Choose 'No' to only display results in the console.")

if (useOutputFolder) {
    outputFolder = QuPathDialogs.promptForDirectory("Select Output Folder", null)
    
    if (outputFolder == null) {
        def continueWithoutOutput = Dialogs.showYesNoDialog("Continue?", 
            "No output folder selected. Continue with console output only?")
        if (!continueWithoutOutput) {
            println("Analysis cancelled by user.")
            return
        }
    } else {
        println("\nOutput folder selected: ${outputFolder.getAbsolutePath()}")
        
        def timestamp = new SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date())
        def imageName = getCurrentImageData()?.getServer()?.getMetadata()?.getName() ?: "Unknown_Image"
        imageName = imageName.replaceAll("[^a-zA-Z0-9._-]", "_")
        
        outputFolder = new File(outputFolder, "ColocAnalysis_${imageName}_${timestamp}")
        outputFolder.mkdirs()
        
        println("Results will be saved to: ${outputFolder.getAbsolutePath()}")
    }
}

// Cellpose segmentation

def pixelSizeInput = Dialogs.showInputDialog("Segmentation Parameters", "Pixel size (μm):", "0.5")
def diameterInput = Dialogs.showInputDialog("Segmentation Parameters", "Expected cell diameter (pixels, 0 for auto):", "15")
def modelInput = Dialogs.showInputDialog("Segmentation Parameters", "Cellpose model (cyto, cyto2, cyto3, nuclei, or custom path):", "cyto3")

def saveOverlays = false
def saveSingleCellAnalysis = false
if (outputFolder != null) {
    saveOverlays = Dialogs.showYesNoDialog("Save Options", 
        "Save segmentation overlay images?\n\n" +
        "This will create visual overlays showing:\n" +
        "• Individual channel segmentations\n" +
        "• Combination overlays\n" +
        "• Colocalization heat maps\n\n" +
        "Note: This may take additional time.")
    
    saveSingleCellAnalysis = Dialogs.showYesNoDialog("Single-Cell Analysis", 
        "Perform detailed single-cell colocalization analysis?\n\n" +
        "This will generate:\n" +
        "• Individual cell colocalization profiles\n" +
        "• Cell-by-cell statistics and heat maps\n" +
        "• Detailed per-cell CSV export\n" +
        "• Single-cell distribution plots\n\n" +
        "⚠️  Warning: This may take significantly longer for large datasets.\n" +
        "Recommended for datasets with <500 cells.")
}

def pixelSize = Double.parseDouble(pixelSizeInput ?: "0.5")
def diameter = Double.parseDouble(diameterInput ?: "15")
def modelName = modelInput ?: "cyto3"

// Multi-channel segmentation

println("\n=== STARTING MULTI-CHANNEL SEGMENTATION ===")

def pathObjects = getSelectedObjects()
if (pathObjects.isEmpty()) {
    pathObjects = getAnnotationObjects()
    if (pathObjects.isEmpty()) {
        Dialogs.showErrorMessage("Cellpose", "Please select a parent object or create annotations!")
        return
    }
}

def segmentationResults = [:]
def combinationResults = [:]

// Individual channel segmentation
selectedChannelNamesClean.eachWithIndex { channelName, index ->
    println("Segmenting channel: ${channelName}")
    
    def cellpose = Cellpose2D.builder(modelName)
            .pixelSize(pixelSize)
            .channels(channelName)
            .diameter(diameter)
            .measureShape()
            .measureIntensity()
            .classify("${channelName}_positive")
            .build()
    
    try {
        cellpose.detectObjects(imageData, pathObjects)
        def detectedObjects = getDetectionObjects().findAll { 
            it.getPathClass()?.getName()?.contains("${channelName}_positive") 
        }
        segmentationResults[channelName] = detectedObjects
        println("  - Detected ${detectedObjects.size()} ${channelName}-positive objects")
    } catch (Exception e) {
        println("  - Error segmenting ${channelName}: ${e.getMessage()}")
        segmentationResults[channelName] = []
    }
}

// Combination analysis

println("\n=== ANALYZING CHANNEL COMBINATIONS ===")

def checkOverlap(obj1, obj2) {
    def roi1 = obj1.getROI()
    def roi2 = obj2.getROI()
    
    if (roi1 == null || roi2 == null) return false
    
    def intersection = roi1.getGeometry().intersection(roi2.getGeometry())
    def overlapArea = intersection.getArea()
    def minArea = Math.min((double)roi1.getArea(), (double)roi2.getArea())
    
    return (overlapArea / minArea) > 0.5
}

def generateCombinations(list, size) {
    if (size == 1) return list.collect { [it] }
    if (size == list.size()) return [list]
    if (size > list.size()) return []
    
    def result = []
    for (int i = 0; i <= list.size() - size; i++) {
        def head = list[i]
        def tail = list[(i+1)..-1]
        generateCombinations(tail, size - 1).each { combo ->
            result.add([head] + combo)
        }
    }
    return result
}

def channelCombinations = []
for (int size = 2; size <= selectedChannelNamesClean.size(); size++) {
    channelCombinations.addAll(generateCombinations(selectedChannelNamesClean, size))
}

channelCombinations.each { combo ->
    def comboName = combo.join("+")
    def comboObjects = []
    
    if (combo.size() == 2) {
        // Double positive
        def objects1 = segmentationResults[combo[0]] ?: []
        def objects2 = segmentationResults[combo[1]] ?: []
        
        objects1.each { obj1 ->
            objects2.each { obj2 ->
                if (checkOverlap(obj1, obj2)) {
                    obj1.setPathClass(getPathClass("${comboName}_positive"))
                    if (!comboObjects.contains(obj1)) {
                        comboObjects.add(obj1)
                    }
                }
            }
        }
    } else if (combo.size() == 3) {
        // Triple positive - find objects that overlap with ALL three channels
        def objects1 = segmentationResults[combo[0]] ?: []
        def objects2 = segmentationResults[combo[1]] ?: []
        def objects3 = segmentationResults[combo[2]] ?: []
        
        objects1.each { obj1 ->
            def overlapsWithCh2 = false
            def overlapsWithCh3 = false
            
            // Check if obj1 overlaps with any object from channel 2
            objects2.each { obj2 ->
                if (checkOverlap(obj1, obj2)) {
                    overlapsWithCh2 = true
                }
            }
            
            // Check if obj1 overlaps with any object from channel 3
            objects3.each { obj3 ->
                if (checkOverlap(obj1, obj3)) {
                    overlapsWithCh3 = true
                }
            }
            
            // If obj1 overlaps with both other channels, it's triple positive
            if (overlapsWithCh2 && overlapsWithCh3) {
                obj1.setPathClass(getPathClass("${comboName}_positive"))
                if (!comboObjects.contains(obj1)) {
                    comboObjects.add(obj1)
                }
            }
        }
    } else if (combo.size() > 3) {
        // Multi-positive (4+ channels) - generalized approach
        def referenceObjects = segmentationResults[combo[0]] ?: []
        
        referenceObjects.each { refObj ->
            def overlapsWithAll = true
            
            // Check if reference object overlaps with all other channels
            for (int i = 1; i < combo.size(); i++) {
                def channelObjects = segmentationResults[combo[i]] ?: []
                def overlapsWithThisChannel = false
                
                channelObjects.each { obj ->
                    if (checkOverlap(refObj, obj)) {
                        overlapsWithThisChannel = true
                    }
                }
                
                if (!overlapsWithThisChannel) {
                    overlapsWithAll = false
                    break
                }
            }
            
            if (overlapsWithAll) {
                refObj.setPathClass(getPathClass("${comboName}_positive"))
                if (!comboObjects.contains(refObj)) {
                    comboObjects.add(refObj)
                }
            }
        }
    }
    
    combinationResults[comboName] = comboObjects
    println("${comboName}: ${comboObjects.size()} objects")
}

// Colocalization analysis

println("\n=== COLOCALIZATION ANALYSIS ===")

def performColocalizationAnalysis(objects, ch1Index, ch2Index, analysisName) {
    if (objects.isEmpty()) {
        println("No objects for ${analysisName} analysis")
        return null
    }
    
    def pearsonValues = []
    def m1Values = []
    def m2Values = []
    def intensityRatios = []
    def ch1Intensities = []
    def ch2Intensities = []
    def icqValues = []
    def vanSteenselValues = []
    
    def ch1Background = 100
    def ch2Background = 100
    
    objects.each { object ->
        try {
            def roi = object.getROI()
            def request = RegionRequest.createInstance(getCurrentImageData().getServer().getPath(), 1.0, roi)
            def pathImage = IJTools.convertToImagePlus(getCurrentImageData().getServer(), request)
            def imp = pathImage.getImage()
            
            def ch1Image = imp.getStack().getProcessor(ch1Index).convertToFloatProcessor()
            def ch2Image = imp.getStack().getProcessor(ch2Index).convertToFloatProcessor()
            def ch1Pixels = ch1Image.getPixels()
            def ch2Pixels = ch2Image.getPixels()
            
            def mask = createObjectMask(pathImage, object, "cell").getPixels()
            
            def ch1Values = []
            def ch2Values = []
            for (int i = 0; i < ch1Pixels.size(); i++) {
                if (mask[i]) {
                    ch1Values.add(ch1Pixels[i])
                    ch2Values.add(ch2Pixels[i])
                }
            }
            
            if (ch1Values.size() == 0) return
            
            def ch1Mean = ch1Values.sum() / ch1Values.size()
            def ch2Mean = ch2Values.sum() / ch2Values.size()
            
            // Standard Pearson correlation
            def pearsonNum = 0.0
            def ch1SumSq = 0.0
            def ch2SumSq = 0.0
            
            ch1Values.eachWithIndex { val1, i ->
                def val2 = ch2Values[i]
                def diff1 = val1 - ch1Mean
                def diff2 = val2 - ch2Mean
                pearsonNum += diff1 * diff2
                ch1SumSq += diff1 * diff1
                ch2SumSq += diff2 * diff2
            }
            
            def pearson = 0
            if (ch1SumSq > 0 && ch2SumSq > 0) {
                pearson = pearsonNum / Math.sqrt(ch1SumSq * ch2SumSq)
            }
            
            // Li's Intensity Correlation Quotient (ICQ)
            def icq = calculateICQ(ch1Values, ch2Values, ch1Mean, ch2Mean)
            
            // Van Steensel's Cross-correlation
            def vanSteensel = calculateVanSteenselCorrelation(ch1Values, ch2Values, ch1Mean, ch2Mean)
            
            // Standard Manders coefficients
            def m1Num = 0.0, m1Den = 0.0, m2Num = 0.0, m2Den = 0.0
            ch1Values.eachWithIndex { val1, i ->
                def val2 = ch2Values[i]
                if (val2 > ch2Background) m1Num += Math.max(val1 - ch1Background, 0)
                if (val1 > ch1Background) m2Num += Math.max(val2 - ch2Background, 0)
                m1Den += Math.max(val1 - ch1Background, 0)
                m2Den += Math.max(val2 - ch2Background, 0)
            }
            
            def m1 = m1Den > 0 ? m1Num / m1Den : 0.0
            def m2 = m2Den > 0 ? m2Num / m2Den : 0.0
            
            pearsonValues.add(pearson)
            m1Values.add(m1)
            m2Values.add(m2)
            intensityRatios.add(ch2Mean / ch1Mean)
            ch1Intensities.add(ch1Mean)
            ch2Intensities.add(ch2Mean)
            icqValues.add(icq)
            vanSteenselValues.add(vanSteensel)
            
            object.getMeasurementList().putMeasurement("Pearson_${analysisName}", pearson)
            object.getMeasurementList().putMeasurement("M1_${analysisName}", m1)
            object.getMeasurementList().putMeasurement("M2_${analysisName}", m2)
            object.getMeasurementList().putMeasurement("ICQ_${analysisName}", icq)
            object.getMeasurementList().putMeasurement("VanSteensel_${analysisName}", vanSteensel)
            object.getMeasurementList().putMeasurement("IntensityRatio_${analysisName}", ch2Mean / ch1Mean)
            
        } catch (Exception e) {
            println("Error analyzing object in ${analysisName}: ${e.getMessage()}")
        }
    }
    
    return [
        pearson: pearsonValues,
        m1: m1Values,
        m2: m2Values,
        intensityRatio: intensityRatios,
        ch1Intensity: ch1Intensities,
        ch2Intensity: ch2Intensities,
        icq: icqValues,
        vanSteensel: vanSteenselValues,
        objectCount: objects.size()
    ]
}

def calculateICQ(ch1Values, ch2Values, ch1Mean, ch2Mean) {
    // Li's Intensity Correlation Quotient
    // ICQ = (number of pixels with (I1-mean1)*(I2-mean2) > 0) / total pixels - 0.5
    // Range: -0.5 (complete segregation) to +0.5 (complete colocalization)
    
    def positiveProduct = 0
    def totalPixels = ch1Values.size()
    
    ch1Values.eachWithIndex { val1, i ->
        def val2 = ch2Values[i]
        if ((val1 - ch1Mean) * (val2 - ch2Mean) > 0) {
            positiveProduct++
        }
    }
    
    return ((double)positiveProduct / (double)totalPixels) - 0.5
}

def calculateVanSteenselCorrelation(ch1Values, ch2Values, ch1Mean, ch2Mean) {
    // Van Steensel's Cross-correlation with shift analysis
    // This implementation calculates correlation at shift=0 and compares with shifts
    // Returns the correlation coefficient that accounts for systematic offsets
    
    def maxShift = Math.min(5, (int)(ch1Values.size() / 10)) // Limit shift to reasonable range
    def correlations = []
    
    // Calculate correlations for different shifts
    for (int shift = -maxShift; shift <= maxShift; shift++) {
        def correlation = 0
        def validPixels = 0
        def ch1SumSq = 0
        def ch2SumSq = 0
        def crossSum = 0
        
        for (int i = 0; i < ch1Values.size(); i++) {
            def j = i + shift
            if (j >= 0 && j < ch2Values.size()) {
                def ch1Diff = ch1Values[i] - ch1Mean
                def ch2Diff = ch2Values[j] - ch2Mean
                
                crossSum += ch1Diff * ch2Diff
                ch1SumSq += ch1Diff * ch1Diff
                ch2SumSq += ch2Diff * ch2Diff
                validPixels++
            }
        }
        
        if (ch1SumSq > 0 && ch2SumSq > 0 && validPixels > 0) {
            correlation = crossSum / Math.sqrt(ch1SumSq * ch2SumSq)
            correlations.add([shift: shift, correlation: correlation])
        }
    }
    
    // Return the maximum correlation found (indicating best alignment)
    if (correlations.isEmpty()) return 0
    
    def maxCorrelation = correlations.max { it.correlation }
    return maxCorrelation.correlation
}

def colocalizationResults = [:]

// Pairwise colocalization analysis
for (int i = 0; i < selectedChannelNumbers.size(); i++) {
    for (int j = i + 1; j < selectedChannelNumbers.size(); j++) {
        def ch1Name = selectedChannelNamesClean[i]
        def ch2Name = selectedChannelNamesClean[j]
        def ch1Num = selectedChannelNumbers[i]
        def ch2Num = selectedChannelNumbers[j]
        def analysisName = "${ch1Name}_vs_${ch2Name}"
        
        def analysisObjects = combinationResults["${ch1Name}+${ch2Name}"] ?: 
                            (segmentationResults[ch1Name] + segmentationResults[ch2Name]).unique()
        
        def result = performColocalizationAnalysis(analysisObjects, ch1Num, ch2Num, analysisName)
        if (result != null) {
            colocalizationResults[analysisName] = result
        }
    }
}

// Triple-positive colocalization analysis (if 3+ channels selected)
if (selectedChannelNumbers.size() >= 3) {
    println("\n=== TRIPLE-POSITIVE COLOCALIZATION ANALYSIS ===")
    
    // For each triple combination, perform pairwise analysis within the triple-positive population
    def tripleCombinations = generateCombinations(selectedChannelNamesClean, 3)
    
    tripleCombinations.each { triple ->
        def tripleComboName = triple.join("+")
        def tripleObjects = combinationResults[tripleComboName] ?: []
        
        if (!tripleObjects.isEmpty()) {
            println("Analyzing triple-positive population: ${tripleComboName} (${tripleObjects.size()} cells)")
            
            // Perform pairwise analysis within triple-positive cells
            for (int i = 0; i < 3; i++) {
                for (int j = i + 1; j < 3; j++) {
                    def ch1Name = triple[i]
                    def ch2Name = triple[j]
                    def ch1Num = selectedChannelNumbers[selectedChannelNamesClean.indexOf(ch1Name)]
                    def ch2Num = selectedChannelNumbers[selectedChannelNamesClean.indexOf(ch2Name)]
                    def analysisName = "${tripleComboName}_${ch1Name}_vs_${ch2Name}"
                    
                    def result = performColocalizationAnalysis(tripleObjects, ch1Num, ch2Num, analysisName)
                    if (result != null) {
                        colocalizationResults[analysisName] = result
                        println("  ${ch1Name} vs ${ch2Name} in triple-positive cells: ${result.objectCount} analyzed")
                    }
                }
            }
        } else {
            println("No triple-positive cells found for: ${tripleComboName}")
        }
    }
}

// Statistical reporting
println("\n" + "=" * 80)
println("COMPREHENSIVE COLOCALIZATION ANALYSIS REPORT")
println("=" * 80)

println("\nSEGMENTATION SUMMARY:")
println("-" * 40)
selectedChannelNamesClean.each { channelName ->
    def count = segmentationResults[channelName]?.size() ?: 0
    println("${channelName}: ${count} objects detected")
}

println("\nCOMBINATION ANALYSIS:")
println("-" * 40)
combinationResults.each { comboName, objects ->
    println("${comboName}: ${objects.size()} co-positive objects")
}

println("\nCOLOCALIZATION STATISTICS:")
println("-" * 40)

colocalizationResults.each { analysisName, results ->
    println("\n${analysisName.toUpperCase()} (n=${results.objectCount}):")
    
    def pearsonStats = calculateSummaryStats(results.pearson)
    def m1Stats = calculateSummaryStats(results.m1)
    def m2Stats = calculateSummaryStats(results.m2)
    def ratioStats = calculateSummaryStats(results.intensityRatio)
    
    println("  Pearson Correlation:")
    println("    Mean ± SD: ${String.format('%.3f', pearsonStats.mean)} ± ${String.format('%.3f', pearsonStats.std)}")
    println("    Median [Range]: ${String.format('%.3f', pearsonStats.median)} [${String.format('%.3f', pearsonStats.min)} - ${String.format('%.3f', pearsonStats.max)}]")
    
    println("  Manders M1 (Ch1 in Ch2 areas):")
    println("    Mean ± SD: ${String.format('%.3f', m1Stats.mean)} ± ${String.format('%.3f', m1Stats.std)}")
    
    println("  Manders M2 (Ch2 in Ch1 areas):")
    println("    Mean ± SD: ${String.format('%.3f', m2Stats.mean)} ± ${String.format('%.3f', m2Stats.std)}")
    
    println("  Intensity Ratio (Ch2/Ch1):")
    println("    Mean ± SD: ${String.format('%.2f', ratioStats.mean)} ± ${String.format('%.2f', ratioStats.std)}")
    
    def colocalizationLevel = ""
    if (pearsonStats.mean > 0.7) colocalizationLevel = "Strong positive"
    else if (pearsonStats.mean > 0.3) colocalizationLevel = "Moderate positive"
    else if (pearsonStats.mean > -0.3) colocalizationLevel = "Weak/No"
    else colocalizationLevel = "Negative"
    
    println("  Interpretation: ${colocalizationLevel} colocalization")
}

println("\n" + "=" * 80)
println("CSV EXPORT DATA:")
println("=" * 80)

colocalizationResults.each { analysisName, results ->
    println("\n${analysisName}_SUMMARY:")
    println("Metric,Mean,Std,Median,Min,Max,Count")
    
    def pearsonStats = calculateSummaryStats(results.pearson)
    def m1Stats = calculateSummaryStats(results.m1)
    def m2Stats = calculateSummaryStats(results.m2)
    def icqStats = calculateSummaryStats(results.icq)
    def vanSteenselStats = calculateSummaryStats(results.vanSteensel)
    
    println("Pearson,${String.format('%.4f', (double)pearsonStats.mean)},${String.format('%.4f', (double)pearsonStats.std)},${String.format('%.4f', (double)pearsonStats.median)},${String.format('%.4f', (double)pearsonStats.min)},${String.format('%.4f', (double)pearsonStats.max)},${results.objectCount}")
    println("M1,${String.format('%.4f', (double)m1Stats.mean)},${String.format('%.4f', (double)m1Stats.std)},${String.format('%.4f', (double)m1Stats.median)},${String.format('%.4f', (double)m1Stats.min)},${String.format('%.4f', (double)m1Stats.max)},${results.objectCount}")
    println("M2,${String.format('%.4f', (double)m2Stats.mean)},${String.format('%.4f', (double)m2Stats.std)},${String.format('%.4f', (double)m2Stats.median)},${String.format('%.4f', (double)m2Stats.min)},${String.format('%.4f', (double)m2Stats.max)},${results.objectCount}")
    println("ICQ,${String.format('%.4f', (double)icqStats.mean)},${String.format('%.4f', (double)icqStats.std)},${String.format('%.4f', (double)icqStats.median)},${String.format('%.4f', (double)icqStats.min)},${String.format('%.4f', (double)icqStats.max)},${results.objectCount}")
    println("VanSteensel,${String.format('%.4f', (double)vanSteenselStats.mean)},${String.format('%.4f', (double)vanSteenselStats.std)},${String.format('%.4f', (double)vanSteenselStats.median)},${String.format('%.4f', (double)vanSteenselStats.min)},${String.format('%.4f', (double)vanSteenselStats.max)},${results.objectCount}")
}

println("\nAnalysis complete! Results have been added to object measurements and can be exported via QuPath's measurement export function.")

// Save results
if (outputFolder != null) {
    println("\n" + "=" * 80)
    println("SAVING RESULTS TO OUTPUT FOLDER")
    println("=" * 80)
    
    try {
        saveHTMLReport(outputFolder, selectedChannelNamesClean, segmentationResults, 
                      combinationResults, colocalizationResults, pixelSize, diameter, modelName)
        
        saveCSVResults(outputFolder, colocalizationResults)
        
        saveObjectMeasurements(outputFolder, colocalizationResults)
        
        saveSummaryStatistics(outputFolder, selectedChannelNamesClean, segmentationResults, 
                             combinationResults, colocalizationResults)
        
        if (saveSingleCellAnalysis) {
            performSingleCellAnalysis(outputFolder, selectedChannelNamesClean, segmentationResults, 
                                    colocalizationResults, selectedChannelNumbers)
        }
        
        if (saveOverlays) {
            println("Creating overlay images...")
            println("Overlay creation feature coming soon!")
        }
        
        println("\nAll results saved successfully to:")
        println(outputFolder.getAbsolutePath())
        
    } catch (Exception e) {
        println("Error saving results: ${e.getMessage()}")
        e.printStackTrace()
    }
}
