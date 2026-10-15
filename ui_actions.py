class EventHandler( object ):
    def doMouseDown( self, evt ):
        evt.Skip()
        return self
    def doMouseUp( self, evt ):
        evt.Skip()
        return self
    def doMouseMove( self, evt ):
        evt.Skip()
        return self
    def doMouseDrag( self, evt ):
        evt.Skip()
        return self
    
    def doKeyDown( self, evt ):
        evt.Skip()
        return self
    def doKeyUp( self, evt ):
        evt.Skip()
        return self
    
    def doContextMenu( self, evt ):
        evt.Skip()
        return self
    
    def doDraw( self ):
        return self


def forward_events( from_event_handler_or_canvas, to_event_handler ):
    assert from_event_handler_or_canvas is not to_event_handler
    
    canvas = from_event_handler_or_canvas
    event_handler = to_event_handler
    
    if not hasattr( canvas, '_doMouseDown' ): canvas._doMouseDown = canvas.doMouseDown
    if not hasattr( canvas, '_doMouseUp' ): canvas._doMouseUp = canvas.doMouseUp
    if not hasattr( canvas, '_doMouseMove' ): canvas._doMouseMove = canvas.doMouseMove
    if not hasattr( canvas, '_doMouseDrag' ): canvas._doMouseDrag = canvas.doMouseDrag
    if not hasattr( canvas, '_doKeyDown' ): canvas._doKeyDown = canvas.doKeyDown
    if not hasattr( canvas, '_doKeyUp' ): canvas._doKeyUp = canvas.doKeyUp
    if not hasattr( canvas, '_doContextMenu' ): canvas._doContextMenu = canvas.doContextMenu
    
    canvas.doMouseDown = lambda *args, **kwargs: event_handler.doMouseDown( *args, **kwargs )
    canvas.doMouseUp = lambda *args, **kwargs: event_handler.doMouseUp( *args, **kwargs )
    canvas.doMouseMove = lambda *args, **kwargs: event_handler.doMouseMove( *args, **kwargs )
    canvas.doMouseDrag = lambda *args, **kwargs: event_handler.doMouseDrag( *args, **kwargs )
    canvas.doKeyDown = lambda *args, **kwargs: event_handler.doKeyDown( *args, **kwargs )
    canvas.doKeyUp = lambda *args, **kwargs: event_handler.doKeyUp( *args, **kwargs )
    canvas.doContextMenu = lambda *args, **kwargs: event_handler.doContextMenu( *args, **kwargs )
    
    assert not hasattr( event_handler, '_event_parent' )
    event_handler._event_parent = canvas
    assert not hasattr( canvas, '_event_child' )
    canvas._event_child = event_handler

def unforward_events( event_handler ):
    if not hasattr( event_handler, '_event_parent' ):
        print('unforward_events(): event_handler has no parent!')
        return
    
    event_parent = event_handler._event_parent
    event_parent.doMouseDown = event_parent._doMouseDown
    event_parent.doMouseUp = event_parent._doMouseUp
    event_parent.doMouseMove = event_parent._doMouseMove
    event_parent.doMouseDrag = event_parent._doMouseDrag
    event_parent.doKeyDown = event_parent._doKeyDown
    event_parent.doKeyUp = event_parent._doKeyUp
    event_parent.doContextMenu = event_parent._doContextMenu
    del event_parent._event_child
    del event_handler._event_parent

def event_parent( event_handler ):
    if not hasattr( event_handler, '_event_parent' ):
        print('event_parent(): event_handler has no parent!')
        return None
    else:
        return event_handler._event_parent

def event_child( event_handler ):
    if not hasattr( event_handler, '_event_child' ):
        ## This method shouldn't print anything; it is used to check whether an object is
        ## forwarding its events.
        #print 'event_child(): event_handler has no child!'
        return None
    else:
        return event_handler._event_child

class DebugEventHandler( EventHandler ):
    def doMouseDown( self, evt ):
        print(evt)
        return self
    def doMouseUp( self, evt ):
        print(evt)
        return self
    def doMouseMove( self, evt ):
        print(evt)
        return self
    def doMouseDrag( self, evt ):
        print(evt)
        return self
    
    def doKeyDown( self, evt ):
        print(evt)
        return self
    def doKeyUp( self, evt ):
        print(evt)
        return self
    
    def doContextMenu( self, evt ):
        print(evt)
        return evt
    
    def doDraw( self ):
        print('drawing')
        
        from OpenGL.GL import glClearColor, glClear, GL_COLOR_BUFFER_BIT, GL_DEPTH_BUFFER_BIT
        
        glClearColor( .5, 0, .5 )
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT )
        return self

class StrokeHandler( EventHandler ):
    def __init__( self, list_to_manage, point_added_cb, complete_cb, event2world, debugging = False ):
        
        self.stroke = list_to_manage
        ## So that the user can pass None to list_to_manage and
        ## rely on the parameter to complete_cb
        if self.stroke is None: self.stroke = []
        
        self.point_added_callback = point_added_cb
        self.complete_callback = complete_cb
        
        self.event2world = event2world
        
        self.debugging = debugging
    
    def doMouseDown( self, evt ):
        del self.stroke[:]
        self.stroke.append( self.event2world( evt.GetPosition() ) )
        if self.debugging: print(self.stroke[-1])
        if self.point_added_callback is not None: self.point_added_callback( self.stroke[-1] )
        return self
    
    def doMouseDrag( self, evt ):
        pt = self.event2world( evt.GetPosition() )
        
        if pt != self.stroke[-1]:
            self.stroke.append( pt )
            if self.debugging: print(self.stroke[-1])
            if self.point_added_callback is not None: self.point_added_callback( self.stroke[-1] )
        
        return self
    
    def doMouseUp( self, evt ):
        pt = self.event2world( evt.GetPosition() )
        
        if pt != self.stroke[-1]:
            self.stroke.append( pt )
            if self.debugging: print(self.stroke[-1])
            if self.point_added_callback is not None: self.point_added_callback( self.stroke[-1] )
        
        if self.complete_callback is not None: self.complete_callback( self.stroke )
        return None

class RotateHandler( EventHandler ):
    def __init__( self, camera, moved_cb, event2screen, debugging = False ):
        
        self.camera = camera
        
        self.moved_callback = moved_cb
        
        self.event2screen = event2screen
        
        self.last_pos = None
        
        self.debugging = debugging
    
    def doMouseDown( self, evt ):
        pos = self.event2screen( evt.GetPosition() )
        
        if self.debugging: print(pos)
        
        self.last_pos = pos
        return self
    
    def doMouseDrag( self, evt ):
        pos = self.event2screen( evt.GetPosition() )
        if pos != self.last_pos:
            self.camera.Rotate2( *[ -( pos[i] - self.last_pos[i] ) for i in range(2) ] )
            #print "camera.eye:", self.camera.eye
            if self.moved_callback is not None: self.moved_callback()
        
        if self.debugging: print(pos)
        
        self.last_pos = pos
        return self
    
    def doMouseUp( self, evt ):
        return self.doMouseDrag( evt )
