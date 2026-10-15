## Everything I need from shapes and manips should go here
class Primitive( object ):
    pass

## All annotations should go here
class Annotation( object ):
    pass

class Graph( object ):
    def __init__( self, anns ):
        self.anns = anns
    
    def node2children( self, node ):
        '''
        Returns a list of Primitives.
        '''
        return []
    def node2connections( self, node ):
        '''
        Returns a list of connections.
        '''
        return []
    def node2connection_children( self, node ):
        '''
        Returns a list of ( connection, children ) tuples.
        '''
        return []
