from wx_annotate_collection import *

import cuboid.shape

class DrawCuboidMode( ModeProxy ):
    def __init__( self, canvas ):
        ModeProxy.__init__( self, canvas )
    
    def doMouseDown( self, evt ):
        ## forward events to a spine drawer.
        ## TODO Q: As soon as finished drawing, switch to ArrowMode?
        forward_events( self, CuboidAnnotationCreateMode.DrawCurveCreator( self ) )
        self.doMouseDown( evt )
    
    def doDraw( self ):
        
        def draw_label():
            self.prepare2D()
            draw_mode_label_text( 'Draw a curve' )
        
        callback = draw_label
        
        child = event_child( self )
        if child is not None:
            callback = lambda: child.doDraw()
        
        self.DrawScene( callback = callback )

class CuboidAnnotationCreateMode( ModeProxy ):
    SAMPLES_PER_UNIT_DISTANCE = 20
    
    class DrawCurveCreator( StrokeDrawer ):
        def __init__( self, canvas ):
            StrokeDrawer.__init__( self, canvas, [] )
            ## Any lookups which might pass through to our proxied class must be done after
            ## calling the superclass constructor.
        
        def doPointAdded( self ):
            self.Refresh()
        
        def doStrokeFinished( self ):
            if len( self.stroke ) >= 3:
                
                resample_N = max( 3, int(ceil( CuboidAnnotationCreateMode.SAMPLES_PER_UNIT_DISTANCE * self.GetZoom() * line_strip_length( self.stroke ) )) )
                
                if self.GetMouseState().ShiftDown():
                    self.stroke = asarray(
                        lerps(
                            self.stroke[0],
                            self.stroke[-1],
                            linspace( 0, 1, resample_N )
                            )
                        )
                else:
                    self.stroke = smooth_open_stroke(
                        self.stroke,
                        resample_N = resample_N
                        )
                
                shape = cuboid.shape.CuboidShape.FromCurveAndRadius( self.stroke, GetSweepShapeCurrentRadius( self ) )
                
                shape.update_3d()
                self.SaveStateForUndo()
                self.annotations.anns.append( shape )
                self.annotations.Beautify()
                self.GetTopLevelParent().Refresh()
            
            #self.shape = None
            self.Refresh()
            
            ## Switch to arrow mode.
            #self.SetModeByID( self.GetModeIDForShortName( "Select" ) )
            self.SwitchToolbarTools( self.GetModeIDForShortName( "Select" ) )
            self.mode.selection.add( shape )
        
        def doSetGLLineStyle( self ):
            glColor4ub( *( Color[ 'spine' ] + (200,) ) )
            self.glLineWidthWorld( self.doGetLineWidthWorld() )
    
    class SpineRedrawer( StrokeDrawer ):
        def __init__( self, canvas, shape ):
            self._shape = shape
            
            StrokeDrawer.__init__( self, canvas, [], auto_unforward_events = False )
            ## Any lookups which might pass through to our proxied class must be done after
            ## calling the superclass constructor.
        
        def get_shape( self ):
            ## If the user performs an undo or redo, self.shape will no longer be in the list
            ## of annotations.  Find the shape in the list of annotations with the same uuid.
            return GetCurrentShapeState( self.annotations.anns, self._shape )
        
        shape = property( get_shape )
        
        def doPointAdded( self ):
            self.Refresh()
        
        def doStrokeFinished( self ):
            if len( self.stroke ) >= 2:
                
                spine_N = len( self.shape.spine2d )
                
                stroke = smooth_open_stroke(
                    self.stroke,
                    resample_N = max( 3, int(ceil( SweepShapeAnnotationCreateMode.SAMPLES_PER_UNIT_DISTANCE * line_strip_length( self.stroke ) )) )
                    )
                
                self.SaveStateForUndo()
                import curve_oversketching
                #self.shape.set_spine2d( stroke )
                old_ts = self.shape.get_all_particular_ts()
                new_ts = list( old_ts )
                ## If the oversketch curve and spine2d are not sampled uniformly, t's along the
                ## spine will be broken by any subsequent resampling;
                ## SweepShape's resample() and smooth_spine() methods are sensitive to this.
                new_spine = curve_oversketching.curve_from_base_curve_and_oversketch_curve(
                    self.shape.spine2d, stroke,
                    scale = self.GetZoom(),
                    ts = new_ts
                    )
                assert len( new_ts ) == len( old_ts )
                self.shape.set_spine2d( new_spine, dict( list(zip( old_ts, new_ts )) ) )
                self.shape.resample( max(
                        spine_N,
                        max( 3, int(ceil( SweepShapeAnnotationCreateMode.SAMPLES_PER_UNIT_DISTANCE * line_strip_length( new_spine ) )) )
                        ) )
                ## Smooth again, for good measure.
                ## UPDATE: Actually, smoothing on curves is pretty heavy; too heavy for oversketching.
                #self.shape.smooth_spine2d()
                self.shape.update_3d()
                self.annotations.Beautify()
                self.GetTopLevelParent().Refresh()
                #print "successfully oversketched"
            
            else:
                ## A click exits the mode.
                unforward_events( self )
            
            self.Refresh()
        
        def doSetGLLineStyle( self ):
            glColor4ub( *( Color[ 'spine' ] + (200,) ) )
            self.glLineWidthWorld( self.doGetLineWidthWorld() )
        
        def doDraw( self ):
            if self.shape is None:
                unforward_events( self )
                self.Refresh()
                return
            
            self.prepare2D()
            draw_mode_label_text( 'Oversketch' )
            StrokeDrawer.doDraw( self )
    
    class CircularCrossSectionTilter( StrokeDrawer ):
        def __init__( self, canvas, shape, t = None ):
            ## The implementation of this class assumes that the spine exists.
            assert len( shape.spine2d ) >= 2
            
            StrokeDrawer.__init__( self, canvas, [] )
            self.SaveStateForUndo()
            self.shape = shape
            self.t = t
            if self.t is not None:
                assert self.t in self.shape.tilt_radians
                assert self.t >= 0.
                assert self.t <= 1.
        
        def doPointAdded( self ):
            ## Is this the first point?  If so, pick the closest spine point.
            if self.t is None:
                assert len( self.stroke ) == 1
                distsqr, pt, edge_index = closest_distsqr_and_point_and_edge_index_on_line_strip_to_point(
                    self.shape.spine2d, self.stroke[0]
                    )
                ## Make sure the closest edge doesn't already have a tilt.
                ## TODO Q: Should I assign one regardless, and leave it up to the shape object
                ##         to tessellate itself sufficiently?
                if distsqr < EdgeDistanceThreshold( self )**2 and edge_index not in [ self.shape.edge_index_at_t( t ) for t in self.shape.tilt_radians.keys() ]:
                    self.t = annotations_sweep.edge_index2t( edge_index, self.shape.spine2d )
                else:
                    unforward_events( self )
                    return
            
            ## We need to update tilt radians.
            self.shape.tilt_radians[ self.t ] = self.shape.compute_tilt_radians_from_point_at_t( self.t, self.stroke[-1] )
            
            if kRealTimeHacks:
                ## Find the shape with the same uuid and change its tilt radians.
                beauty = [
                    shape for shape in self.annotations.anns_beautified.of_type( SweepShape )
                    if shape.uuid == self.shape.uuid
                    ]
                ## There should be one shape with the same uuid.
                assert len( beauty ) == 1
                beauty = beauty[0]
                beauty.tilt_radians[ self.t ] = self.shape.tilt_radians[ self.t ]
                ## If no-pinch is enabled, temporarily turn it off and generate the mesh.
                no_pinch = annotations_sweep.kNoPinch
                if no_pinch:
                    annotations_sweep.kNoPinch = False
                    beauty.get_mesh()
                    annotations_sweep.kNoPinch = True
                self.GetTopLevelParent().Refresh()
            
            self.Refresh()
        
        def doStrokeFinished( self ):
            ## If we were updating the mesh in real-time, we temporarily turned of no-pinch.
            ## Now regenerate the mesh with no-pinch enabled.
            if kRealTimeHacks and annotations_sweep.kNoPinch:
                self.shape.gen_mesh()
            
            self.shape.update_3d()
            self.annotations.Beautify()
            self.GetTopLevelParent().Refresh()
            
            self.shape = None
        
        def doDraw( self ):
            ## We override doDraw() so that we don't draw the stroke.  Instead, draw the circle.
            ## Actually, adding the tilted circle to self.shape.tilted_circles means that the
            ## regular shape drawing code draws the cross section.
            self.DrawCuboidShapeHandles( [ self.shape ] )
    
    class CircularCrossSectionRadiusManip( StrokeDrawer ):
        def __init__( self, canvas, shape, t = None ):
            ## The implementation of this class assumes that the spine exists.
            assert len( shape.spine2d ) >= 2
            
            StrokeDrawer.__init__( self, canvas, [] )
            self.SaveStateForUndo()
            self.shape = shape
            self.t = t
            if self.t is not None:
                assert self.t in self.shape.radii
                assert self.t >= 0.
                assert self.t <= 1.
        
        def doPointAdded( self ):
            ## Is this the first point?  If so, pick the closest spine point.
            if self.t is None:
                assert len( self.stroke ) == 1
                ## Cast a ray perpendicular to the closest point on the shape's boundary?
                offsets = annotations_sweep.boundary_offsets2d( self.shape.spine2d, self.shape.dense_radii )
                spine2d = asarray( self.shape.spine2d )
                distsqr, pt, edge_index = min(
                    closest_distsqr_and_point_and_edge_index_on_line_strip_to_point(
                        spine2d + offsets, self.stroke[0]
                        ),
                    closest_distsqr_and_point_and_edge_index_on_line_strip_to_point(
                        spine2d - offsets, self.stroke[0]
                        )
                    )
                ## Make sure the closest edge doesn't already have a cross section.
                ## TODO Q: Should I assign one regardless, and leave it up to the shape object
                ##         to tessellate itself sufficiently?
                if distsqr < EdgeDistanceThreshold( self )**2 and edge_index not in [ self.shape.edge_index_at_t( t ) for t in self.shape.radii.keys() ]:
                    self.t = annotations_sweep.edge_index2t( edge_index, self.shape.spine2d )
                    ## Initialize the cross section to a list of one item.
                    ## (The list's item will be set properly below.)
                    self.shape.radii[ self.t ] = None
                else:
                    unforward_events( self )
                    return
            
            ## Update the circle's radius.
            center = self.shape.spine2d_at_t( self.t )
            radius = mag( asarray( self.stroke[-1] ) - center )
            self.shape.radii[ self.t ] = radius
            
            self.Refresh()
        
        def doStrokeFinished( self ):
            self.shape.update_3d()
            self.annotations.Beautify()
            self.GetTopLevelParent().Refresh()
            
            self.shape = None
        
        def doDraw( self ):
            ## We override doDraw() so that we don't draw the stroke.  Instead, draw the circle.
            ## Actually, adding the tilted circle to self.shape.tilted_circles means that the
            ## regular shape drawing code draws the cross section.
            self.DrawCuboidShapeHandles( [ self.shape ] )
    
    class Translater( StrokeDrawer ):
        def __init__( self, canvas, shape ):
            ## The implementation of this class assumes that the spine exists.
            assert len( shape.spine2d ) >= 2
            
            StrokeDrawer.__init__( self, canvas, [] )
            self.shape = shape
        
        def doPointAdded( self ):
            if len( self.stroke ) < 2: return
            elif len( self.stroke ) == 2: self.SaveStateForUndo()
            
            delta = asarray( self.stroke[-1] ) - asarray( self.stroke[-2] )
            self.shape.translate_2d( delta )
            
            self.Refresh()
        
        def doStrokeFinished( self ):
            if len( self.stroke ) >= 2:
                ## Actually, I don't need to call this for translate_2d().
                #self.shape.update_3d()
                self.annotations.Beautify()
                self.GetTopLevelParent().Refresh()
            else:
                self.Refresh()
            
            self.shape = None
        
        def doDraw( self ):
            ## Don't draw anything, we're translating the shape, so normal drawing shows all.
            pass
    
    def doMouseDown( self, evt ):
        ## This is crucial for contextual menus to work.  Mouse down gets first dibs on a potential
        ## contextual menu event, but the right mouse button will be down.
        if evt.RightDown():
            print('skipping mouse down')
            evt.Skip()
            return
        
        ## If we pick a handle (circular cross section or curve to overdraw), manipulate/overdraw it.
        ## Otherwise, create a new annotation.
        
        picked_handle = None
        parts = None
        ## If the command key is down, always draw a new spine.
        if evt.CmdDown():
            parts = []
        ## If the shift key is down, always move.
        elif evt.ShiftDown():
            parts = [ CuboidShapePart.BoundaryInterior ]
        picked_handle = PickClosestCuboidShapeHandle(
            self.annotations.anns.of_type( cuboid.shape.CuboidShape ),
            self,
            self.event2world( evt.GetPosition() ),
            parts = parts
            )
        if picked_handle is not None:
            print("XXX TODO CONTINUE: ADAPT ME")
            
            #print "Picked shape, handle:", picked_handle
            if picked_handle[1][0] == SweepShapePart.CrossSectionTiltHandle:
                forward_events( self, self.CircularCrossSectionTilter( self, picked_handle[0], picked_handle[1][1] ) )
                self.doMouseDown( evt )
            elif picked_handle[1][0] == SweepShapePart.RadiusHandle:
                forward_events( self, self.CircularCrossSectionRadiusManip( self, picked_handle[0], picked_handle[1][1] ) )
                self.doMouseDown( evt )
            elif picked_handle[1][0] == SweepShapePart.SymmetryPlaneHandle:
                forward_events( self, self.SymmetryPlaneManip( self, picked_handle[0], picked_handle[1][1] ) )
                self.doMouseDown( evt )
            elif picked_handle[1][0] == SweepShapePart.Spine:
                forward_events( self, self.CircularCrossSectionTilter( self, picked_handle[0] ) )
                self.doMouseDown( evt )
            elif picked_handle[1][0] == SweepShapePart.Boundary:
                forward_events( self, self.CircularCrossSectionRadiusManip( self, picked_handle[0] ) )
                self.doMouseDown( evt )
            elif picked_handle[1][0] in (SweepShapePart.BoundaryInterior,):
                forward_events( self, self.Translater( self, picked_handle[0] ) )
                self.doMouseDown( evt )
        else:
            forward_events( self, self.SpineDrawer( self ) )
            self.doMouseDown( evt )
    
    def doContextMenu( self, evt ):
        ## NOTE: I can't use evt.GetPosition() because it returns bogus data here.
        world_pos = self.event2world( self.GetMousePosition() )
        
        ## If we pick a SweepShape, show the contextual menu.
        picked_handle = PickClosestSweepShapeHandle(
            self.annotations.anns.of_type( SweepShape ),
            self,
            world_pos
            )
        if picked_handle is None:
            evt.Skip()
            return
        
        shape = picked_handle[0]
        
        menu = wx.Menu()
        menu_ids = []
        
        print('TODO CONTINUE')
        
        def redraw_silhouette( evt ):
            forward_events( self, self.SilhouetteRedrawer( self, shape ) )
        
        item = menu.Append( wx.ID_ANY, "Redraw Silhouette" )
        self.Bind( wx.EVT_MENU, redraw_silhouette, id = item.GetId() )
        menu_ids.append( item.GetId() )
        
        def redraw_spine( evt ):
            forward_events( self, self.SpineRedrawer( self, shape ) )
        
        item = menu.Append( wx.ID_ANY, "Redraw Spine" )
        self.Bind( wx.EVT_MENU, redraw_spine, id = item.GetId() )
        menu_ids.append( item.GetId() )
        
        if True: #wx_annotate_shell.kDebugging:
            def smooth():
                self.SaveStateForUndo()
                shape.smooth_spine2d()
                GarbageCollectMetaAnnotations( self.annotations.anns )
                shape.update_3d()
                
                self.annotations.Beautify()
                self.GetTopLevelParent().Refresh()
            
            item = menu.Append( wx.ID_ANY, "Smooth Spine" )
            self.Bind( wx.EVT_MENU, lambda evt: smooth(), id = item.GetId() )
            menu_ids.append( item.GetId() )
            
            
            def resample( factor ):
                self.SaveStateForUndo()
                shape.resample( max( 3, int(round( factor * len( shape.spine2d ) )) ) )
                GarbageCollectMetaAnnotations( self.annotations.anns )
                shape.update_3d()
                
                self.annotations.Beautify()
                self.GetTopLevelParent().Refresh()
            
            item = menu.Append( wx.ID_ANY, "Double Spine Tesselation" )
            self.Bind( wx.EVT_MENU, lambda evt: resample( 2 ), id = item.GetId() )
            menu_ids.append( item.GetId() )
            
            item = menu.Append( wx.ID_ANY, "Halve Spine Tesselation" )
            self.Bind( wx.EVT_MENU, lambda evt: resample( .5 ), id = item.GetId() )
            menu_ids.append( item.GetId() )
        
        ## Did we pick a tilt or radius manipulation handle?  If so, add its items to the contextual menu.
        picked = PickSweepShapeHandle( shape, self, world_pos, parts = (SweepShapePart.CrossSectionTiltHandle,) )
        if picked is not None:
            t = picked[1]
            
            menu.AppendSeparator()
            
            ## Add delete cross section to the menu.
            ## Make delete disabled if we don't want the user to delete the cross section
            ## (start or end cross section or only one left).
            
            def delete_circle( evt, t = t ):
                self.SaveStateForUndo()
                del shape.tilt_radians[ t ]
                GarbageCollectMetaAnnotations( self.annotations.anns )
                shape.update_3d()
                
                self.annotations.Beautify()
                self.GetTopLevelParent().Refresh()
            
            item = menu.Append( wx.ID_ANY, "Delete Tilt Handle" )
            self.Bind( wx.EVT_MENU, delete_circle, id = item.GetId() )
            menu_ids.append( item.GetId() )
            
            ## Don't allow users to delete the start or end cross sections.
            ## UPDATE: Actually, I don't care.
            ## UPDATE 2: Radii falloff to zero or less with only one cross section on either end.
            ## TODO: Fix radii falling-off with only one cross section on either end.
            #if t == 0. or t == 1.: item.Enable( False )
            
            ## Don't allow users to delete the last remaining cross section.
            if len( shape.tilt_radians ) == 1: item.Enable( False )
        
        picked = PickSweepShapeHandle( shape, self, world_pos, parts = (SweepShapePart.RadiusHandle,) )
        if picked is not None:
            t = picked[1]
            
            menu.AppendSeparator()
            
            ## Add delete cross section to the menu.
            ## Make delete disabled if we don't want the user to delete the cross section
            ## (start or end cross section or only one left).
            
            def delete_circle( evt, t = t ):
                self.SaveStateForUndo()
                del shape.radii[ t ]
                GarbageCollectMetaAnnotations( self.annotations.anns )
                shape.update_3d()
                
                self.annotations.Beautify()
                self.GetTopLevelParent().Refresh()
            
            item = menu.Append( wx.ID_ANY, "Delete Scale Handle" )
            self.Bind( wx.EVT_MENU, delete_circle, id = item.GetId() )
            menu_ids.append( item.GetId() )
            
            ## Don't allow users to delete the start or end cross sections.
            ## UPDATE: Actually, I don't care.
            ## UPDATE 2: Radii falloff to zero or less with only one cross section on either end.
            ## TODO: Fix radii falling-off with only one cross section on either end.
            #if t == 0. or t == 1.: item.Enable( False )
            
            ## Don't allow users to delete the last remaining cross section.
            if len( shape.radii ) == 1: item.Enable( False )
        
        ## Did we pick the spine?  If so, allow the user to draw an arbitrary cross section.
        ## If there already is an arbitrary cross section draw, allow the user to delete it.
        ## TODO: Add an arbitrary cross section "pencil" handle for control-clicking on.
        ##       Or this will be handled by proper selection UI for cross sections.
        picked = PickSweepShapeHandle( shape, self, world_pos, parts = (SweepShapePart.Spine,) )
        ## If we don't pick the spine, also check if we picked an existing arbitrary cross section.
        ## TODO: Arbitrary cross sections are not drawn as circles, but they are picked as circles.
        ##       Fix this so that they are picked as arbitrary cross sections.
        if picked is None:
            picked = PickSweepShapeHandle( shape, self, world_pos, parts = (SweepShapePart.CrossSectionInterior, SweepShapePart.CrossSection) )
            ## Make sure this picked cross section is really an arbitrary one.
            if picked is not None and picked[1] not in shape.cross_sections: picked = None
        if picked is not None:
            t = picked[1]
            
            menu.AppendSeparator()
            
            ## Add draw arbitrary cross section to the menu.
            
            def draw_arbitrary_cross_section( evt, t = t ):
                ## Push the new UI mode.
                forward_events( self, self.ArbitraryCrossSection( self, shape, t, immediate_click = False ) )
                
                '''
                ## For now, just insert a dummy cross section.
                self.SaveStateForUndo()
                #shape.cross_sections[ t ] = annotations_sweep.unit_circle
                ## A cardioid.
                r = 1./3.
                shape.cross_sections[ t ] = asarray( [ (
                        2*r * ( cos( theta ) - .5 * cos( 2*theta ) ),
                        2*r * ( sin( theta ) - .5 * sin( 2*theta ) )
                    ) for theta in linspace( 0, 2*pi, kCrossSectionSamples )[:-1] ] )
                del r
                pass
                
                GarbageCollectMetaAnnotations( self.annotations.anns )
                self.annotations.Beautify()
                self.GetTopLevelParent().Refresh()
                '''
            
            if t in shape.cross_sections:
                txt = "Re-Draw Arbitrary Cross Section"
            else:
                txt = "Draw Arbitrary Cross Section"
            item = menu.Append( wx.ID_ANY, txt )
            self.Bind( wx.EVT_MENU, draw_arbitrary_cross_section, id = item.GetId() )
            menu_ids.append( item.GetId() )
            
            
            ## Add smooth and delete cross-section to the menu.
            if t in shape.cross_sections:
                ## Add smooth cross section to the menu.
                def smooth_arbitrary_cross_section( evt, t = t ):
                    self.SaveStateForUndo()
                    shape.cross_sections[ t ] = annotations_sweep.normalize_arbitrary_cross_section(
                        smooth_closed_stroke(
                            shape.cross_sections[ t ],
                            resample_N = len( shape.cross_sections[ t ] )
                            )
                        )
                    GarbageCollectMetaAnnotations( self.annotations.anns )
                    shape.update_3d()
                    
                    self.annotations.Beautify()
                    self.GetTopLevelParent().Refresh()
                
                item = menu.Append( wx.ID_ANY, "Smooth Arbitrary Cross Section" )
                self.Bind( wx.EVT_MENU, smooth_arbitrary_cross_section, id = item.GetId() )
                menu_ids.append( item.GetId() )
                
                ## Add delete cross section to the menu.
                def delete_arbitrary_cross_section( evt, t = t ):
                    self.SaveStateForUndo()
                    del shape.cross_sections[ t ]
                    GarbageCollectMetaAnnotations( self.annotations.anns )
                    shape.update_3d()
                    
                    self.annotations.Beautify()
                    self.GetTopLevelParent().Refresh()
                
                item = menu.Append( wx.ID_ANY, "Delete Arbitrary Cross Section" )
                self.Bind( wx.EVT_MENU, delete_arbitrary_cross_section, id = item.GetId() )
                menu_ids.append( item.GetId() )
        
        self.PopupMenu( menu )
        menu.Destroy()
        ## I'm not sure where I got this idea of unbinding these events.
        ## I guess it's just good cleanup, so we don't leave a dangling bound callback to a menu
        ## event that can no longer be called, and also to free the id's we reserved so we
        ## don't run out of id space, since every time this is called we generate new ids.
        for id in menu_ids: self.Unbind( wx.EVT_MENU, id = id )
    
    def doKeyUp( self, evt ):
        ## Some operations change depending on modifier state.
        self.Refresh()
        evt.Skip()
    
    def doKeyDown( self, evt ):
        #print evt.GetKeyCode()
        keycode = evt.GetKeyCode()
        
        ## Some operations change depending on modifier state.
        self.Refresh()
        
        if keycode in kDeleteKeys:
            
            ## evt.GetPosition() returns garbage.
            world_pos = self.event2world( self.GetMousePosition() )
            
            ## First check if the mouse is over an annotation.
            ## If the mouse is not over any annotations, check if we're over a handle.
            ann = None
            
            ## TODO: Convenience function for this annotation picking.
            #ann = self.PickAnnotationOfType( screen_pos, annotations_ellipsoid.AnnotationData )
            t_anns = PickCuboidShapes(
                self.annotations.anns.of_type( cuboid.shape.CuboidShape ),
                world_pos
                )
            if len( t_anns ) > 0:
                ann = min( t_anns )[1]
            
            if ann is not None:
                self.SaveStateForUndo()
                print('deleting an annotation')
                RemoveAnnotation( self.annotations.anns, ann )
                
                self.annotations.Beautify()
                self.GetTopLevelParent().Refresh()
        
        else:
            evt.Skip()
    
    def doMouseMove( self, evt ):
        ## Refresh so our status text changes.
        self.Refresh()
    
    def doDraw( self ):
        def draw_child():
            child = event_child( self )
            if child is not None:
                child.doDraw()
                return True
        
        def draw_callback():
            if draw_child(): return
            
            self.prepare2D()
            
            text = ''
            
            picked_handle = None
            parts = None
            ## If the command key is down, always draw a new spine.
            if self.GetMouseState().CmdDown():
                parts = []
            ## If the shift key is down, always move.
            elif self.GetMouseState().ShiftDown():
                parts = [ CuboidShapePart.BoundaryInterior ]
            picked_handle = PickClosestCuboidShapeHandle(
                self.annotations.anns.of_type( cuboid.shape.CuboidShape ),
                self,
                self.event2world( self.GetMousePosition() ),
                parts = parts
                )
            if picked_handle is not None:
                print("XXX TODO CONTINUE: ADAPT ME")
                
                draw_handles = [ picked_handle[0] ]
                ## XXX TODO: This is ugly.
                #if wx_annotate_shell.kPaperFigureMode: draw_handles = None
                self.DrawCuboidShapeHandles( draw_handles )
                
                if picked_handle[1][0] == SweepShapePart.CrossSectionTiltHandle:
                    text = 'Change tilt'
                elif picked_handle[1][0] == SweepShapePart.RadiusHandle:
                    ## This text is too long
                    #text = 'Change circular cross section radius'
                    text = 'Change thickness'
                elif picked_handle[1][0] == SweepShapePart.SymmetryPlaneHandle:
                    text = 'Change symmetry plane'
                elif picked_handle[1][0] == SweepShapePart.Spine:
                    text = 'Tilt'
                elif picked_handle[1][0] == SweepShapePart.Boundary:
                    text = 'Change thickness'
                elif picked_handle[1][0] in (SweepShapePart.BoundaryInterior,):
                    text = 'Move'
            else:
                text = 'Draw a spine'
            draw_mode_label_text( text )
        
        self.DrawScene( callback = draw_callback )
#===================
