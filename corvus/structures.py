from abc import ABCMeta, abstractmethod

# Debug: FDV
import sys

# Debug: FDV
import pprint
pp_debug = pprint.PrettyPrinter(indent=4)

class Workflow(object):
    def __init__(self, target, desc=''):

# Modified by FDV
# We have already ensured that things are in the right format, so not much
# testing
#       if isinstance(target, list) and target and isinstance(target[0], basestring):
#           self.target = target
#       elif isinstance(target, basestring):
#           self.target = [target]
#       else:
#           raise TypeError('Target should be token or list of tokens')
# The input formatter generates a list of lists. Here we only use the
# first element to initialize the target.
        self.target = target[0]
        self.sequence = []
        self.cost = 0
        self.desc = desc
        self.userinput = self.target
    
    def __str__(self):
        lines = []
        if len(self.desc) > 0:
            lines.append(self.desc)
        for i, e in enumerate(self.sequence):
            lines.append(str(i+1) + ') ' + str(e))
        lines.append('Required User Input: ' + str(self.userinput))
        return '\n'.join(lines)

    def addExchange(self, exchange):
        self.addExchangeAt(len(self.sequence), exchange) 

    def addExchangeAt(self, index, exchange):
        self.sequence.insert(index, exchange)
        self.updateCost()
        self.updateRequiredInput()
        
    def addExchangeList(self, list):
        self.addExchangeListAt(len(self.sequence), list)
        
    def addExchangeListAt(self, i, list):
        self.sequence[i:i] = list
        self.updateCost()
        self.updateRequiredInput()

    def getCost(self):
        return self.cost

    def getRequiredInput(self):
        return self.userinput

    def updateCost(self):
        self.cost = 0
        for xc in self.sequence:
            self.cost = self.cost + xc.cost
    
    def updateRequiredInput(self):    
        inp = set()
        for s in reversed(self.sequence):
            inp.update(set(s.input))
            inp.difference_update(set(s.output))
        self.userinput = list(inp)

class Exchange(object):
    def __init__(self, handler, input, output, cost=1, desc=''):
        self.handler = handler
        self.input = input
        self.output = output
        self.cost = cost
        self.desc = desc

    def __str__(self):
        lines = [self.desc]
        lines.append('   input: ' + str(self.input))
        lines.append('  output: ' + str(self.output))
        return '\n'.join(lines)

    def go(self, config, system):
        returnedOutput = {}
        for token in self.output:
            returnedOutput[token] = None
# Debug:FDV
#lines below uncommented thru sys.exit, Krsna
        #print self
        #print self.handler
        #print self.handler.exchange
        #print type(self)
        #print type(self.handler)
        #print type(self.handler.exchange)
        #sys.exit()
        self.handler.exchange(config, system, returnedOutput)
#       print('\nsystem\n')
#       pp_debug.pprint(system)
        system.update(returnedOutput)

class Update(object):
    def __init__(self, token, newToken=None, newValue=None):
        isStr = lambda x: isinstance(x, str)
        isList = lambda x: isinstance(x, list)
        isStrList = lambda L: isList(L) and L and all(map(isStr, L))
        if isStr(token):
            # One token updated
            self.output = [token]
            self.desc = 'Update ' + token
            if isStr(newToken) and not newValue:
                self.dynamic = True
                self.input = [newToken]
                self.desc += ' with value of ' + newToken
            elif newValue is not None and not newToken:
                self.dynamic = False
                self.input = [newValue]
                self.desc += ' with value = ' + str(newValue)
            else:
                msg = 'Error creating Update: use either newToken or newValue'
                raise Exception(msg)
        elif isStrList(token):
            # Multiple tokens updated
            self.output = token
            self.desc = 'Update ' + str(token)
            if isStrList(newToken) and len(newToken) == len(token) and not newValue:
                self.dynamic = True
                self.input = newToken
                self.desc += ' with values of ' + str(newToken)
            elif isList(newValue) and len(newValue) == len(token) and not newToken:
                self.dynamic = False
                self.input = newValue
                self.desc += ' with values = ' + str(newValue)
            else:
                msg = 'Error creating Update: use either newToken or newValue'
                raise Exception(msg)
            self.desc += ' respectively'
        else: 
            msg = 'Error creating Update: need token or list of tokens to update'
            raise Exception(msg)
        self.cost = 0

    def __str__(self):
        return self.desc

    def go(self, config, system):
        for i, token in enumerate(self.output):
            if self.dynamic:
                system[token] = system[self.input[i]]
            else:
                system[token] = self.input[i]

class Loop(object):
    def __init__(self, exchange, token, gridToken=None, gridValues=None):
        isStr = lambda x: isinstance(x, str)
        isList = lambda x: isinstance(x, list)
        isStrList = lambda L: isList(L) and L and all(map(isStr, L))
        if isStr(token):
            # Loop over one parameter
            self.param = [token]
            self.desc = 'Loop over ' + token
            if isStr(gridToken) and not gridValues:
                self.dynamic = True
                self.input = [gridToken]
                self.desc += ' with dynamic grid ' + gridToken
            elif isList(gridValues) and not gridToken:
                self.dynamic = False
                self.input = [gridValues]
                self.desc += ' with predefined grid ' + str(gridValues)
            else:
                msg = 'Error creating Update: use either gridToken or gridValues'
                raise Exception(msg)
        elif isStrList(token):
            # Loop over multiple paramenters
            self.param = token
            self.desc = 'Loop over ' + str(token)
            if isStrList(gridToken) and len(gridToken) == len(token) and not gridValues:
                self.dynamic = True
                self.input = gridToken
                self.desc += ' with dynamic grids from ' + str(gridToken)
            elif isList(gridValues) and len(gridValues) == len(token) and not gridToken:
                self.dynamic = False
                self.input = gridValues
                self.desc += ' with predefined grids ' + str(gridValues)
            else:
                msg = 'Error creating Update: use either gridToken or gridValues'
                raise Exception(msg)
            self.desc += ' respectively'
        else:
            msg = 'Error creating Loop: need token or list of tokens to loop over'
            raise Exception(msg)
        self.exchange = exchange
        self.cost = exchange.cost
        self.output = [o + '-grid' for o in exchange.output]

    def __str__(self):
        lines = [self.desc]
        lines.append('   input: ' + str(self.exchange.input))
        lines.append('  output: ' + str(self.output))
        return '\n'.join(lines)

    def go(self, config, system):
        baseIndex = str(config['xcIndex'])
        # Set up loop table of results
        table = {}
        for i,p in enumerate(self.param):
            if self.dynamic:
                table[p + '-grid'] = system[self.input[i]]
            else:
                table[p + '-grid'] = self.input[i]
        # Make sure if multiple grids need to update, they all have the same length
        gridlengths = [len(table[g]) for g in table]
        if gridlengths.count(gridlengths[0]) == len(gridlengths):
            gridlength = gridlengths[0]
        else:
            msg = 'Loop Error: all input grids should be same length'
            msg += str(table)
            raise Exception(msg)
        # Initialize table with space for output grids 
        for token in self.exchange.output:
            table[token + '-grid'] = []
        # Loop through, Update grids, run Exchange
        for i in range(gridlength):
            for p in self.param:
                Update(p, newValue=table[p + '-grid'][i]).go(config, system)
            config['xcIndex'] = baseIndex + '.' + str(i+1)
            self.exchange.go(config, system)
            for token in self.exchange.output:
                table[token + '-grid'].append(system[token])
        system.update(table)                 

class Handler(metaclass=ABCMeta):
    
    @classmethod
    def exchange(self, config, input, output):
        self.prep(config)
# NOTE FDV:
# These changes are to make structure.py compatible with JJK's changes in the
# feff.py handler. He implemented inserting the generateInput andi
# translateOutput methods inside the run method, as we had discussed to allow
# for more generic handlers that don't used file based input output.
# Here I try to reproduce what he did in his on structures.py
#       inputFiles = self.generateInput(config, input, output)
#       self.run(config, inputFiles)
        self.run(config, input, output)
#       self.translateOutput(config, input, output)
        self.cleanup(config)

    @abstractmethod
    def prep(config):
        pass

#   @abstractmethod
#   def generateInput(config, input, output):
#       pass

#   @abstractmethod
#   def run(config, files):
#       pass
    @abstractmethod
    def run(config, input, output):
        pass

#   @abstractmethod
#   def translateOutput(config, input, output):
#       pass

    @abstractmethod
    def cleanup(config):
        pass

    @abstractmethod
    def canProduce(output):
        # Return True if Handler can calculate all output(s)
        # output could be single string or list of strings
        pass

    @abstractmethod
    def requiredInputFor(output):
        # Return list of tokens (strings)
        # output could be single string or list of strings
        pass

    @abstractmethod
    def costOf(output):
        # Return integer indicating costliness of calculation
        # output could be single string or list of strings
        pass

    @abstractmethod
    def sequenceFor(output,inp=None):
        # Return list of Workflow sequence(s) to produce output
        # output could be single string or list of strings
        pass
