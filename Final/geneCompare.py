# Group Members: Fabrice Kurmann, Ronnie Voskoboynik

import statistics
import sys

########################################################################
# CSVReader Reader 
########################################################################
class CSVReader:
    '''A simple CSV file reader that yields each line of a CSV file.'''
    def __init__ (self, fname=None):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is None:
            return sys.stdin
        else:
            return open(self.fname)
        
    def readCSV (self):
        '''Read a CSV file, return it line by line.'''
        lineArray = []
        
        with self.doOpen() as fileH:
            lineArray = []
            # for line in fileH:
            line = fileH.readline()
            while line:
            
                lineArray = line.split(',')
                yield lineArray

                line = fileH.readline()

########################################################################
# GeneStats 
########################################################################
class GeneStats:
    '''Class that performs all information compiling and statistical analysis on 
    gene information from an input file.'''
    def __init__(self, reader):
        self.csvReader = reader

    def sumExpression(self, dictionary):
        '''Fill a dictionary with the total expression of all genes.'''
        expressionDict = dictionary
        # Read in data row by row, skipping first row
        firstRow = 1
        for list in self.csvReader.readCSV():
            if (firstRow == 1):
                firstRow = 0
                continue
            
            # Get the name of the gene
            gene = list[1]
            # Get the solitary expression
            totalSolitary = 0
            for i in range(2,8):
                totalSolitary += int(list[i])
            
            # Get the chimeric expression
            totalChimeric = 0
            for i in range(8,len(list)):
                totalChimeric += int(list[i])

            # Either add to a gene's expression list or add a kv pair to the expression dictionary
            if expressionDict.get(gene):
                expressionDict[gene][0].append(totalSolitary)
                expressionDict[gene][1].append(totalSolitary)
            else:
                expressionDict[gene] = ([totalSolitary], [totalChimeric])
        return expressionDict
    
    def getStatistics(self, expressionDict, dictionary):
        '''Fill a dictionary with the average expression of all genes.'''
        # Find averages of expressions across all genes
        averageDict = dictionary
        for key, value in expressionDict.items():
            averageSolitary, averageChimeric = 0, 0
            # Sum up all values for the gene's solitary and chimeric features
            for i in range(len(value[0])):
                averageSolitary += value[0][i]
                averageChimeric += value[1][i]
            # Compute averages
            averageSolitary = averageSolitary / (6 * len(value[0]))
            averageChimeric = averageChimeric / (6 * len(value[0]))

            # Compute standard deviation of solitary gene values
            # For genes that appear once, no standard deviation can be found, use average * 0.25
            if (len(value[0]) == 1):
                stdSolitary = averageSolitary / 4
            else:
                stdSolitary = statistics.stdev(value[0])

            # Add to average dictionary
            averageDict[key] = (averageSolitary, averageChimeric, stdSolitary)
        return averageDict

########################################################################
# Main 
########################################################################
def main(inFile = None, maxStandardDeviations=1000000000, minStandardDeviations=3, includeZeros=False, numResults = 100):
    '''First create a reader object to read the input file, then create a GeneStats object.
    Then run all the functions necessary to find genes with the greatest variance
    between solitary and chimeric examples.
    
    Usage: 
    Optional arguments:
    inFile: the file to read from, must be CSV
    maxStandardDeviations: the upper limit for number of standard deviations the chimeric genes may be from the solitary
    minStandardDeviations: the lower limit for number of standard deviations the chimeric genes may be from the solitary
    includeZeros: boolean value for whether the returnList includes all genes that are not expressdd
    numResults: the upper limit for number of results that will be returned, results are returned in decreasing tScore order, so smaller values may be cut
'''
    Reader = CSVReader(inFile) 
    GeneInformation = GeneStats(Reader)
    
    # In this execution we are starting with a blank dictionary for each run
    inputDict = {}

    # Format {gene: ([list of solitary expression totals], [list of chimeric expression totals])}
    expressionDict = GeneInformation.sumExpression(inputDict)
    # print(expressionDict)
    
    # Format {gene: (average expression solitary, average expression chimeric)}
    averageDict = GeneInformation.getStatistics(expressionDict, inputDict)
    # print(averageDict)

    # Create a list that will be sorted before its contents are returned
    # Format [(name, averageSolitary, averageChimeric, tScore),...]
    returnList = []
    
    # Output all gene information where solitary and chimeric fall within the specified number of standard deviations
    for key, value in averageDict.items():
        averageSolitary = value[0]
        averageChimeric = value[1]
        standardDeviation = value[2]

        # Show those genes falling under the specified amount of standard deviations
        tScore = 0
        if standardDeviation == 0:
            pass
        else:
            tScore = abs(averageChimeric - averageSolitary) / standardDeviation
        if tScore <= maxStandardDeviations:
            # Only include genes who's solitary value are zero if the user specifies
            if (includeZeros == False and averageSolitary == 0):
                pass
            # Only include genes with a meaninful standard deviation
            elif (tScore < minStandardDeviations):
                pass
            else:
                returnList.append((key, averageSolitary, averageChimeric, tScore))
            
    # Sort the return list based on the highest t-scores
    returnList = sorted(returnList, key=lambda x: -1*x[3])

    # Print results to stdout
    print ('Returning genes withing the given variance threshold, sorted by t-score, descending')
    print ('Format: Gene Name     Solitary Gene Average Expression      Chimeric Gene Average Expression     T-score')
    index = 0
    for item in returnList:
        index += 1
        print ('{:15s} : {:10.2f} {:37.2f} {:34.2f}'.format(item[0], item[1], item[2], item[3]))
        if (index == numResults):
            break
    
if __name__ == "__main__":
    main(inFile='Data.csv', includeZeros=False, maxStandardDeviations=10)
