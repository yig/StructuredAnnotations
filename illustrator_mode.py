from wx_annotate_collection import *

class ArrowMode( ModeProxy ):
    def __init__( self, canvas ):
        ModeProxy.__init__( self, canvas )
        
        self.selection = set()
    
    def FixupSelection( self ):
        '''
        Because undo/redo changes the objects in the list of annotations, the selection can point
        to objects in the undo/redo stack.  Fix them up to point to objects in the current set.
        '''
        
        self.selection = set( [
            GetCurrentShapeState( self.annotations.anns, shape )
            for shape in self.selection
            ] )
        ## If undo/redo included the creation or destruction of a shape, then GetCurrentShapeState()
        ## will return None.  Remove all None's from the selection.
        try: self.selection.remove( None )
        except KeyError: pass
    
    def PickShape( self, world_pos, world_dir = (0,0,-1) ):
        ## Skip locked shapes.
        return self.PickShapeFromSequence(
            [ shape for shape in self.annotations.anns if shape not in self.locked_shapes ],
            world_pos, world_dir
            )
    
    def PickShapeFromSequence( self, shapes, world_pos, world_dir = (0,0,-1), shapes_are_selected = False ):
        if shapes_are_selected:
            ## TODO Q: Should I just make parts None, so that anything is picked?
            ##         The danger with an explicit list is that I forget to update it.
            parts = [SweepShapePart.CrossSectionTiltHandle,SweepShapePart.RadiusHandle,SweepShapePart.SymmetryPlaneHandle,SweepShapePart.CapRadiusHandle,SweepShapePart.CrossSection,SweepShapePart.Spine,SweepShapePart.Boundary,SweepShapePart.BoundaryInterior]
            kind2world_radius = lambda kind: PointDistanceThreshold( self )
        else:
            parts = [SweepShapePart.BoundaryInterior]
            kind2world_radius = lambda kind: 0.
        
        for shape in shapes:
            if isinstance( shape, SweepShape ):
                handle = PickSweepShapeHandle( shape, self, world_pos, world_dir, parts = parts )
                if handle is not None:
                    assert handle[0] in parts
                    return shape
            
            elif isinstance( shape, annotations_ellipsoid.AnnotationData ):
                ## TODO: pick only the circular cross section
                handle = wx_annotate_sphere.PickAnnotationHandle(
                    shape,
                    world_pos,
                    kind2world_radius
                    )
                if handle is not None:
                    return shape
            
            else:
                ## Is it even a shape?  It could be any kind of annotation.
                #raise NotImplementedError, "No Manip for the shape type"
                pass
        
        return None
    
    def PickHandleFromShapes( self, shapes, world_pos, world_dir = (0,0,-1) ):
        for shape in shapes:
            if isinstance( shape, SweepShape ):
                parts = [SweepShapePart.CrossSectionTiltHandle,SweepShapePart.RadiusHandle,SweepShapePart.SymmetryPlaneHandle,SweepShapePart.CapRadiusHandle,SweepShapePart.CrossSection,SweepShapePart.Spine,SweepShapePart.Boundary,SweepShapePart.BoundaryInterior]
                handle = PickSweepShapeHandle( shape, self, world_pos, world_dir, parts = parts )
                if handle is not None:
                    assert handle[0] in parts
                    return shape, handle
            
            elif isinstance( shape, annotations_ellipsoid.AnnotationData ):
                ## TODO: pick only the circular cross section
                handle = wx_annotate_sphere.PickAnnotationHandle(
                    shape,
                    world_pos,
                    lambda kind: PointDistanceThreshold( self )
                    )
                if handle is not None:
                    return shape, handle
            
            else:
                ## Is it even a shape?  It could be any kind of annotation.
                #raise NotImplementedError, "No Manip for the shape type"
                pass
        
        return None
    
    def TextFromShapeHandle( self, shape_handle ):
        shape, handle = shape_handle
        if isinstance( shape, SweepShape ):
            part2text = {
                SweepShapePart.CrossSectionTiltHandle: 'Change tilt',
                SweepShapePart.RadiusHandle: 'Change thickness',
                SweepShapePart.CapRadiusHandle: 'Change cap size',
                SweepShapePart.SymmetryPlaneHandle: 'Change symmetry plane',
                SweepShapePart.Spine: 'Tilt',
                SweepShapePart.Boundary: 'Change thickness',
                SweepShapePart.BoundaryInterior: 'Move',
                SweepShapePart.CrossSection: 'Edit cross section'
                }
            ## Someone may call this with a handle for which there is no action;
            ## in that case, print nothing.
            return part2text.get( handle[0], '' )
            
        elif isinstance( shape, annotations_ellipsoid.AnnotationData ):
            kind2text = {
                annotations_ellipsoid.AnnotationData.WholeEllipseHandle: 'Move',
                annotations_ellipsoid.AnnotationData.Axis0Handle: 'Rotate/Scale',
                annotations_ellipsoid.AnnotationData.Axis1Handle: 'Change circle radius',
                annotations_ellipsoid.AnnotationData.CrossSectionHandle: 'Tilt circular cross section'
                }
            return kind2text[ handle.kind ]
    
    def ModeFromShapeHandle( self, shape_handle ):
        shape, handle = shape_handle
        if isinstance( shape, SweepShape ):
            if handle[0] == SweepShapePart.CrossSectionTiltHandle:
                return SweepShapeAnnotationCreateMode.CircularCrossSectionTilter( self, shape, handle[1] )
            elif handle[0] == SweepShapePart.RadiusHandle:
                return SweepShapeAnnotationCreateMode.CircularCrossSectionRadiusManip( self, shape, handle[1] )
            elif handle[0] == SweepShapePart.CapRadiusHandle:
                return SweepShapeAnnotationCreateMode.CapRadiusManip( self, shape, handle[1] )
            elif handle[0] == SweepShapePart.SymmetryPlaneHandle:
                return SweepShapeAnnotationCreateMode.SymmetryPlaneManip( self, shape, handle[1] )
            elif handle[0] == SweepShapePart.Spine:
                return SweepShapeAnnotationCreateMode.CircularCrossSectionTilter( self, shape )
                #return SweepShapeSpineDeformer( self, shape )
            elif handle[0] == SweepShapePart.Boundary:
                return SweepShapeAnnotationCreateMode.CircularCrossSectionRadiusManip( self, shape )
            elif handle[0] in (SweepShapePart.BoundaryInterior,):
                return SweepShapeAnnotationCreateMode.Translater( self, shape )
            elif handle[0] in (SweepShapePart.CrossSection,):
                return SweepShapeAnnotationCreateMode.ArbitraryCrossSection( self, shape, handle[1] )
            else:
                raise RuntimeError("Unknown handle type!")
        
        elif isinstance( shape, annotations_ellipsoid.AnnotationData ):
            return EllipseAnnotationCreateMode.HandleDragger( self, handle )
    
    def doMouseDown( self, evt ):
        ## TODO: If undo/redo is called, we have no callback to change the selection.
        self.FixupSelection()
        
        ## I need mouse down so we can manipulate handles.
        #def doMouseUp( self, evt ):
        #if not a click: return
        
        ## This is crucial for contextual menus to work.  Mouse down gets first dibs on a potential
        ## contextual menu event, but the right mouse button will be down.
        if evt.RightDown():
            print('skipping mouse down')
            evt.Skip()
            return
        
        ## If shift is down, we're adding to or (if also option key is down) removing from
        ## the selection.
        if evt.ShiftDown():
            if evt.AltDown():
                try:
                    self.selection.remove( self.PickShape( self.event2world( evt.GetPosition() ) ) )
                except KeyError: pass
            else:
                pick_shape = self.PickShape( self.event2world( evt.GetPosition() ) )
                self.selection.add( pick_shape )
                self.annotations.MoveToFront( pick_shape )
        
        ## Else if command key is down, we're translating.
        elif len( self.selection ) > 0 and evt.CmdDown():
            forward_events( self, ArrowMode.TranslateMultiple( self, self.selection ) )
            self.doMouseDown( evt )
        
        else:
            ## Else if there's a handle under the mouse, we're manipulating that handle.
            shape_handle = self.PickHandleFromShapes( self.selection, self.event2world( evt.GetPosition() ) )
            if shape_handle is not None:
                forward_events( self, self.ModeFromShapeHandle( shape_handle ) )
                self.doMouseDown( evt )
            
            ## Otherwise, we're changing the selection to the object under the mouse.
            else:
                pick_shape = self.PickShape( self.event2world( evt.GetPosition() ) )
                ## Move-to-front when clicking.
                if pick_shape is not None:
                    self.annotations.MoveToFront( pick_shape )
                #print 'PickShape:', pick_shape
                self.selection = set([ pick_shape ])
                try:
                    self.selection.remove( None )
                except KeyError: pass
        
        self.GetTopLevelParent().Refresh()
    
    def doContextMenu( self, evt ):
        #print self.selection
        
        menu = wx.Menu()
        menu_ids = []
        
        ## NOTE: I can't use evt.GetPosition() because it returns bogus data here.
        world_pos = self.event2world( self.GetMousePosition() )
        
        shape = self.PickShapeFromSequence( self.selection, world_pos, shapes_are_selected = True )
        if shape is not None:
            if isinstance( shape, SweepShape ):
                self._AddSweepShapeMenuItems( menu, menu_ids, shape, world_pos )
                menu.AppendSeparator()
            elif isinstance( shape, annotations_ellipsoid.AnnotationData ):
                self._AddEllipseMenuItems( menu, menu_ids, shape, world_pos )
                menu.AppendSeparator()
            else:
                print("No contextual menu items for shape type")
        
        if shape is not None:
            def lock_shape( evt ):
                self.SaveStateForUndo()
                self.locked_shapes.add( shape )
                self.selection.remove( shape )
                self.Refresh()
            item = menu.Append( wx.ID_ANY, "Lock Shape" )
            self.Bind( wx.EVT_MENU, lock_shape, id = item.GetId() )
            menu_ids.append( item.GetId() )
        
        if len( self.selection ) > 0:
            def lock_selected_shapes( evt ):
                self.SaveStateForUndo()
                self.locked_shapes.update( self.selection )
                self.selection.clear()
                self.Refresh()
            item = menu.Append( wx.ID_ANY, "Lock All Selected Shapes" )
            self.Bind( wx.EVT_MENU, lock_selected_shapes, id = item.GetId() )
            menu_ids.append( item.GetId() )
        
        def unlock_all_shapes( evt ):
            self.SaveStateForUndo()
            self.locked_shapes.clear()
            self.Refresh()
        item = menu.Append( wx.ID_ANY, "Unlock All Shapes" )
        self.Bind( wx.EVT_MENU, unlock_all_shapes, id = item.GetId() )
        menu_ids.append( item.GetId() )
        if len( self.locked_shapes ) == 0: item.Enable( False )
        
        
        self.PopupMenu( menu )
        menu.Destroy()
        ## I'm not sure where I got this idea of unbinding these events.
        ## I guess it's just good cleanup, so we don't leave a dangling bound callback to a menu
        ## event that can no longer be called, and also to free the id's we reserved so we
        ## don't run out of id space, since every time this is called we generate new ids.
        for id in menu_ids: self.Unbind( wx.EVT_MENU, id = id )
    
    def _AddSweepShapeMenuItems( self, menu, menu_ids, shape, world_pos ):
        
        ## Gather picking information.
        picked_spine = PickSweepShapeHandle( shape, self, world_pos, parts = (SweepShapePart.Spine,) )
        picked_xsection = PickSweepShapeHandle( shape, self, world_pos, parts = (SweepShapePart.CrossSectionInterior, SweepShapePart.CrossSection) )
        picked_tilt = PickSweepShapeHandle( shape, self, world_pos, parts = (SweepShapePart.CrossSectionTiltHandle,) )
        picked_radius = PickSweepShapeHandle( shape, self, world_pos, parts = (SweepShapePart.RadiusHandle,) )
        picked_symmetry = PickSweepShapeHandle( shape, self, world_pos, parts = (SweepShapePart.SymmetryPlaneHandle,) )
        picked_cap_radius = PickSweepShapeHandle( shape, self, world_pos, parts = (SweepShapePart.CapRadiusHandle,) )
        
        def get_thing_info( picked_thing, thing ):
            if picked_thing is not None:
                return ( True, picked_thing[1] )
            elif picked_spine is not None:
                spine_vertex = shape.vertex_index_at_t( picked_spine[1] )
                collide_ts = [ t_cand for t_cand in thing.keys() if shape.vertex_index_at_t( t_cand ) == spine_vertex ]
                if len( collide_ts ) > 0:
                    return ( True, min( [ ( abs( picked_spine[1] - t ), t ) for t in collide_ts ] )[1] )
                else:
                    return ( False, picked_spine[1] )
            else:
                return None
        
        xsection_info = get_thing_info( picked_xsection, shape.cross_sections )
        tilt_info = get_thing_info( picked_tilt, shape.tilt_radians )
        radius_info = get_thing_info( picked_radius, shape.radii )
        symmetry_info = get_thing_info( picked_symmetry, shape.symmetry_radians )
        ## cap_radius_info is special.
        cap_radius_info = None
        if picked_cap_radius is not None:
            cap_radius_info = [ None, picked_cap_radius[1] ]
        elif picked_spine is not None:
            spine_vertex = shape.vertex_index_at_t( picked_spine[1] )
            if spine_vertex in ( 0, len(shape.spine2d)-1 ):
                cap_radius_info = [ None, 0 if spine_vertex == 0 else 1 ]
        if cap_radius_info is not None:
            cap_radius_info[0] = shape.cap_radii[ cap_radius_info[1] ] is not None
        #print "cap_radius_info:", cap_radius_info
        
        ## First, operations for changing the entire shape.
        def deform_spine( evt ):
            forward_events( self, SweepShapeSpineDeformer( self, shape ) )
        item = menu.Append( wx.ID_ANY, "Deform Spine" )
        self.Bind( wx.EVT_MENU, deform_spine, id = item.GetId() )
        menu_ids.append( item.GetId() )
        
        def redraw_spine( evt ):
            forward_events( self, SweepShapeAnnotationCreateMode.SpineRedrawer( self, shape ) )
        item = menu.Append( wx.ID_ANY, "Over-sketch Spine" )
        self.Bind( wx.EVT_MENU, redraw_spine, id = item.GetId() )
        menu_ids.append( item.GetId() )
        
        def deform_silhouette( evt ):
            forward_events( self, SweepShapeSilhouetteManip( self, shape ) )
        item = menu.Append( wx.ID_ANY, "Deform Silhouette" )
        self.Bind( wx.EVT_MENU, deform_silhouette, id = item.GetId() )
        menu_ids.append( item.GetId() )
        
        def redraw_silhouette( evt ):
            forward_events( self, SweepShapeAnnotationCreateMode.SilhouetteRedrawer( self, shape ) )
        item = menu.Append( wx.ID_ANY, "Redraw Silhouette" )
        ## TODO: Why is this broken? Did it ever work?
        item.Enable( False )
        self.Bind( wx.EVT_MENU, redraw_silhouette, id = item.GetId() )
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
            
            def straighten():
                self.SaveStateForUndo()
                ts = linspace( 0, 1, len( shape.spine2d ) )
                shape.set_spine2d( asarray(
                    lerps( shape.spine2d[0], shape.spine2d[-1], ts )
                    ) )
                GarbageCollectMetaAnnotations( self.annotations.anns )
                shape.update_3d()
                
                self.annotations.Beautify()
                self.GetTopLevelParent().Refresh()
            
            item = menu.Append( wx.ID_ANY, "Straighten Spine" )
            self.Bind( wx.EVT_MENU, lambda evt: straighten(), id = item.GetId() )
            menu_ids.append( item.GetId() )
        
        
        ## Next, control for directly changing radius (this is only useful if we're
        ## ArrowMode2()).
        if isinstance( self, ArrowMode2 ) and picked_radius is not None:
            t = picked_radius[1]
            
            menu.AppendSeparator()
            
            def radius_manip( evt, t = t ):
                forward_events( self, SweepShapeAnnotationCreateMode.CircularCrossSectionRadiusManip( self, shape, t ) )
            item = menu.Append( wx.ID_ANY, "Change Scale" )
            self.Bind( wx.EVT_MENU, radius_manip, id = item.GetId() )
            menu_ids.append( item.GetId() )
        
        
        ## Next, arbitrary cross section controls.
        
        ## TODO: Add an arbitrary cross section "pencil" handle for control-clicking on.
        ##       Or this will be handled by proper selection UI for cross sections.
        if xsection_info is not None:
            t = xsection_info[1]
            
            menu.AppendSeparator()
            
            ## Add draw arbitrary cross section to the menu.
            
            def draw_arbitrary_cross_section( evt ):
                ## Push the new UI mode.
                #forward_events( self, SweepShapeAnnotationCreateMode.ArbitraryCrossSectionDrawer( self, shape, t ) )
                forward_events( self, SweepShapeAnnotationCreateMode.ArbitraryCrossSection( self, shape, t, immediate_click = False ) )
                
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
        
        ## Next, add/remove tilt/radius/symmetry/cap handles.
        if True in [ info is not None for info in [ tilt_info, radius_info, cap_radius_info, symmetry_info ] ]:
            menu.AppendSeparator()
        
        ## tilt handle
        if tilt_info is not None:
            t = tilt_info[1]
            
            if t in shape.tilt_radians:
                ## Add delete cross section to the menu.
                ## Make delete disabled if we don't want the user to delete the cross section
                ## (start or end cross section or only one left).
                
                ## We want the 't' as it is now, not once it is re-defined later.
                ## To capture 't', bind it as a default parameter (from
                ##  http://www.devshed.com/c/a/Python/Nested-Functions-in-Python/1/
                ## ).
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
            else:
                def add_tilt_handle( evt, t = t ):
                    self.SaveStateForUndo()
                    shape.tilt_radians[ t ] = shape.tilt_radians_at_t( t )
                    GarbageCollectMetaAnnotations( self.annotations.anns )
                    shape.update_3d()
                    
                    self.annotations.Beautify()
                    self.GetTopLevelParent().Refresh()
                
                item = menu.Append( wx.ID_ANY, "Add Tilt Handle" )
                self.Bind( wx.EVT_MENU, add_tilt_handle, id = item.GetId() )
                menu_ids.append( item.GetId() )
        
        if radius_info is not None:
            t = radius_info[1]
            
            if t in shape.radii:
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
            else:
                def add_radius_handle( evt, t = t ):
                    self.SaveStateForUndo()
                    shape.radii[ t ] = shape.radius_at_t( t )
                    GarbageCollectMetaAnnotations( self.annotations.anns )
                    shape.update_3d()
                    
                    self.annotations.Beautify()
                    self.GetTopLevelParent().Refresh()
                
                item = menu.Append( wx.ID_ANY, "Add Scale Handle" )
                self.Bind( wx.EVT_MENU, add_radius_handle, id = item.GetId() )
                menu_ids.append( item.GetId() )
        
        if cap_radius_info is not None:
            t = cap_radius_info[1]
            
            ## Add delete cross section to the menu.
            ## Make delete disabled if we don't want the user to delete the cross section
            ## (start or end cross section or only one left).
            
            def delete_cap_radius_handle( evt, t = t ):
                self.SaveStateForUndo()
                shape.cap_radii[ t ] = None
                GarbageCollectMetaAnnotations( self.annotations.anns )
                shape.update_3d()
                
                self.annotations.Beautify()
                self.GetTopLevelParent().Refresh()
            
            item = menu.Append( wx.ID_ANY, "Free Cap Size Handle" )
            self.Bind( wx.EVT_MENU, delete_cap_radius_handle, id = item.GetId() )
            menu_ids.append( item.GetId() )
            
            if not cap_radius_info[0]: item.Enable( False )
        
        if symmetry_info is not None:
            t = symmetry_info[1]
            
            if t in shape.symmetry_radians:
                ## Add "delete symmetry plane handle" to the menu.
                
                def delete_symmetry_plane_handle( evt, t = t ):
                    self.SaveStateForUndo()
                    del shape.symmetry_radians[ t ]
                    GarbageCollectMetaAnnotations( self.annotations.anns )
                    shape.update_3d()
                    
                    self.annotations.Beautify()
                    self.GetTopLevelParent().Refresh()
                
                item = menu.Append( wx.ID_ANY, "Delete Symmetry Plane Handle" )
                self.Bind( wx.EVT_MENU, delete_symmetry_plane_handle, id = item.GetId() )
                menu_ids.append( item.GetId() )
                
                ## Don't allow users to delete the last remaining symmetry radians handle.
                if len( shape.symmetry_radians ) == 1: item.Enable( False )
            else:
                ## Add "add symmetry plane handle", "add radius handle", and "add tilt handle" to the menu.
                
                def add_symmetry_plane_handle( evt, t = t ):
                    self.SaveStateForUndo()
                    shape.symmetry_radians[ t ] = shape.symmetry_radians_at_t( t )
                    GarbageCollectMetaAnnotations( self.annotations.anns )
                    shape.update_3d()
                    
                    self.annotations.Beautify()
                    self.GetTopLevelParent().Refresh()
                
                item = menu.Append( wx.ID_ANY, "Add Symmetry Plane Handle" )
                self.Bind( wx.EVT_MENU, add_symmetry_plane_handle, id = item.GetId() )
                menu_ids.append( item.GetId() )
    
    def _AddEllipseMenuItems( self, menu, menu_ids, shape, world_pos ):
        pass
    
    def doKeyDown( self, evt ):
        ## TODO: If undo/redo is called, we have no callback to change the selection.
        self.FixupSelection()
        
        ## delete the selection
        
        #print evt.GetKeyCode()
        keycode = evt.GetKeyCode()
        
        if keycode in kDeleteKeys:
            
            if len( self.selection ) > 0:
                self.SaveStateForUndo()
                print('deleting selected annotations')
                for shape in self.selection:
                    RemoveAnnotation( self.annotations.anns, shape )
                    try: self.locked_shapes.remove( shape )
                    except KeyError: pass
                self.selection = set()
                
                self.annotations.Beautify()
                self.GetTopLevelParent().Refresh()
        
        else:
            ## If it's a modifier key, it might change what is shown.
            self.Refresh()
            evt.Skip()
    
    def doKeyUp( self, evt ):
        ## Refresh so our status text changes.
        self.Refresh()
    
    def doMouseMove( self, evt ):
        ## Refresh so our status text changes.
        self.Refresh()
    
    def doDraw( self ):
        ## TODO: If undo/redo is called, we have no callback to change the selection.
        self.FixupSelection()
        
        ## Draw shapes, and for shapes in self.selection, draw their handles, too.
        
        def draw_callback():
            ## An attempt to move-to-front when hovering.
            pick_shape = self.PickShape( self.event2world( self.GetMousePosition() ) )
            if pick_shape is not None and ( pick_shape in self.selection or len( self.selection ) == 0 ):
                self.annotations.MoveToFront( pick_shape )
            
            self.DrawSweepShapeHandles( self.selection )
            self.DrawEllipseHandles( self.selection )
            
            ## Draw a picking preview (smart guide).
            if pick_shape is not None and pick_shape not in self.selection:
                if isinstance( pick_shape, SweepShape ):
                    DrawPickedSweepShapeHandle( pick_shape, (SweepShapePart.BoundaryInterior,) )
                elif isinstance( pick_shape, annotations_ellipsoid.AnnotationData ):
                    ## TODO: Smart guides for ellipses.
                    pass
            
            self.prepare2D()
            
            text = ''
            
            if self.GetMouseState().ShiftDown():
                if self.GetMouseState().AltDown():
                    text = 'Remove from selection'
                else:
                    text = 'Add to selection'
            
            ## Else if command key is down, we're translating.
            elif len( self.selection ) > 0 and self.GetMouseState().CmdDown():
                text = 'Translate shape' + ( 's' if len( self.selection ) > 1 else '' )
            
            else:
                ## Else if there's a handle under the mouse, we're manipulating that handle.
                shape_handle = self.PickHandleFromShapes( self.selection, self.event2world( self.GetMousePosition() ) )
                ## Draw Smart Guides.
                if shape_handle is not None:
                    
                    ## Draw smart guides.
                    if isinstance( shape_handle[0], SweepShape ):
                        DrawPickedSweepShapeHandle( shape_handle[0], shape_handle[1] )
                    elif isinstance( shape_handle[0], annotations_ellipsoid.AnnotationData ):
                        ## TODO: Smart guides for ellipses.
                        pass
                    
                    text = self.TextFromShapeHandle( shape_handle )
                
                ## Otherwise, we're changing the selection to the object under the mouse.
                else:
                    text = 'Select shape'
            
            draw_mode_label_text( text )
        
        callback = draw_callback
        
        child = event_child( self )
        if child is not None:
            callback = lambda: child.doDraw()
        
        self.DrawScene( callback = callback )
    
    
    class TranslateMultiple( StrokeDrawer ):
        def __init__( self, canvas, shapes ):
            StrokeDrawer.__init__( self, canvas, [] )
            self.shapes = list(shapes)
        
        def doPointAdded( self ):
            if len( self.stroke ) < 2: return
            elif len( self.stroke ) == 2: self.SaveStateForUndo()
            
            delta = asarray( self.stroke[-1] ) - asarray( self.stroke[-2] )
            for shape in self.shapes: shape.translate( delta )
            
            self.Refresh()
        
        def doStrokeFinished( self ):
            if len( self.stroke ) >= 2:
                ## Actually, I don't need to call this for translate_2d().
                #self.shape.update_3d()
                self.annotations.Beautify()
                self.GetTopLevelParent().Refresh()
            else:
                self.Refresh()
            
            self.shapes = None
        
        def doDraw( self ):
            ## Don't draw anything, we're translating the shape, so normal drawing shows all.
            pass

class ArrowMode2( ArrowMode ):
    def TextFromShapeHandle( self, shape_handle ):
        shape, handle = shape_handle
        if isinstance( shape, SweepShape ):
            part2text = {
                SweepShapePart.CrossSectionTiltHandle: 'Stretch' if self.GetMouseState().AltDown() else 'Change tilt',
                SweepShapePart.RadiusHandle: 'Change scale' if self.GetMouseState().AltDown() else 'Stretch Silhouette',
                #SweepShapePart.RadiusHandle: 'Change thickness',
                #SweepShapePart.RadiusHandle: 'Stretch Silhouette',
                SweepShapePart.CapRadiusHandle: 'Change cap size',
                SweepShapePart.SymmetryPlaneHandle: 'Change symmetry plane',
                SweepShapePart.Spine: 'Tilt' if self.GetMouseState().AltDown() else 'Stretch',
                #SweepShapePart.Spine: 'Tilt',
                #SweepShapePart.Spine: 'Stretch',
                SweepShapePart.Boundary: 'Change scale' if self.GetMouseState().AltDown() else 'Stretch Silhouette',
                #SweepShapePart.Boundary: 'Change thickness',
                #SweepShapePart.Boundary: 'Stretch Silhouette',
                SweepShapePart.BoundaryInterior: 'Move',
                SweepShapePart.CrossSection: 'Edit cross section'
                }
            ## Someone may call this with a handle for which there is no action;
            ## in that case, print nothing.
            return part2text.get( handle[0], '' )
            
        elif isinstance( shape, annotations_ellipsoid.AnnotationData ):
            kind2text = {
                annotations_ellipsoid.AnnotationData.WholeEllipseHandle: 'Move',
                annotations_ellipsoid.AnnotationData.Axis0Handle: 'Rotate/Scale',
                annotations_ellipsoid.AnnotationData.Axis1Handle: 'Change circle radius',
                annotations_ellipsoid.AnnotationData.CrossSectionHandle: 'Tilt circular cross section'
                }
            return kind2text[ handle.kind ]
    
    def ModeFromShapeHandle( self, shape_handle ):
        shape, handle = shape_handle
        if isinstance( shape, SweepShape ):
            if handle[0] == SweepShapePart.CrossSectionTiltHandle:
                if self.GetMouseState().AltDown():
                    return SweepShapeSpineDeformer( self, shape )
                else:
                    return SweepShapeAnnotationCreateMode.CircularCrossSectionTilter( self, shape, handle[1] )
            elif handle[0] == SweepShapePart.RadiusHandle:
                if self.GetMouseState().AltDown():
                    return SweepShapeAnnotationCreateMode.CircularCrossSectionRadiusManip( self, shape, handle[1] )
                else:
                    import silhouette_manip
                    return silhouette_manip.SweepShapeSilhouetteManip( self, shape, handle[1], handle[2] )
            elif handle[0] == SweepShapePart.CapRadiusHandle:
                return SweepShapeAnnotationCreateMode.CapRadiusManip( self, shape, handle[1] )
            elif handle[0] == SweepShapePart.SymmetryPlaneHandle:
                return SweepShapeAnnotationCreateMode.SymmetryPlaneManip( self, shape, handle[1] )
            elif handle[0] == SweepShapePart.Spine:
                if self.GetMouseState().AltDown():
                    return SweepShapeAnnotationCreateMode.CircularCrossSectionTilter( self, shape )
                else:
                    return SweepShapeSpineDeformer( self, shape )
            elif handle[0] == SweepShapePart.Boundary:
                if self.GetMouseState().AltDown():
                    return SweepShapeAnnotationCreateMode.CircularCrossSectionRadiusManip( self, shape )
                else:
                    import silhouette_manip
                    return silhouette_manip.SweepShapeSilhouetteManip( self, shape )
            elif handle[0] in (SweepShapePart.BoundaryInterior,):
                return SweepShapeAnnotationCreateMode.Translater( self, shape )
            elif handle[0] in (SweepShapePart.CrossSection,):
                return SweepShapeAnnotationCreateMode.ArbitraryCrossSection( self, shape, handle[1] )
            else:
                raise RuntimeError("Unknown handle type!")
        
        elif isinstance( shape, annotations_ellipsoid.AnnotationData ):
            return EllipseAnnotationCreateMode.HandleDragger( self, handle )

class DrawSweepShapeMode( ModeProxy ):
    def __init__( self, canvas ):
        ModeProxy.__init__( self, canvas )
    
    def doMouseDown( self, evt ):
        ## forward events to a spine drawer.
        ## TODO Q: As soon as finished drawing, switch to ArrowMode?
        forward_events( self, SweepShapeAnnotationCreateMode.SpineDrawer( self ) )
        self.doMouseDown( evt )
    
    def doDraw( self ):
        
        def draw_label():
            self.prepare2D()
            draw_mode_label_text( 'Draw a spine' )
        
        callback = draw_label
        
        child = event_child( self )
        if child is not None:
            callback = lambda: child.doDraw()
        
        self.DrawScene( callback = callback )


## See "cuboid.ui" for DrawCuboidMode
# class DrawCuboidMode( ModeProxy )

class DrawEllipsoidShapeMode( ModeProxy ):
    def __init__( self, canvas ):
        ModeProxy.__init__( self, canvas )
    
    def doMouseDown( self, evt ):
        ## forward events to ellipsoid drawer.
        ## TODO Q: As soon as finished drawing, switch to ArrowMode?
        forward_events( self, EllipseAnnotationCreateMode.EllipseSilhouetteDrawer( self ) )
        self.doMouseDown( evt )
    
    def doDraw( self ):
        
        def draw_label():
            self.prepare2D()
            draw_mode_label_text( 'Draw an ellipse' )
        
        callback = draw_label
        
        child = event_child( self )
        if child is not None:
            callback = lambda: child.doDraw()
        
        self.DrawScene( callback = callback )


class ZoomMode( ModeProxy ):
    def __init__( self, canvas ):
        ModeProxy.__init__( self, canvas )
        
        self.down_event_loc = None
        
        ## TODO: Change the cursor to a magnifying glass.
        ## TODO: Preserve the selection.
    
    def doMouseDown( self, evt ):
        #print "ZoomMode.doMouseDown()"
        self.down_event_loc = evt.GetPosition()
    
    def doMouseDrag( self, evt ):
        ## So we can draw the rect.
        self.Refresh()
    
    def doMouseUp( self, evt ):
        ## TODO Q: Why is wxWidgets eating double-clicks?
        #print "ZoomMode.doMouseUp()"
        
        ## if click, zoom in or out at the mouse location depending on whether option key is down.
        if self.IsClick( self.down_event_loc, evt.GetPosition() ):
            if evt.AltDown():
                self.ZoomOut( self.event2world( evt.GetPosition() ) )
            else:
                self.ZoomIn( self.event2world( evt.GetPosition() ) )
        
        ## otherwise, zoom so the view contains the rect formed by the mouse-down to mouse-up.
        else:
            self.ZoomToRect( self.event2world( self.down_event_loc ), self.event2world( evt.GetPosition() ) )
        
        self.down_event_loc = None
        self.GetTopLevelParent().Refresh()
    
    def IsClick( self, down_event_pos, up_event_pos ):
        down_event_pos = asarrayf( down_event_pos )
        up_event_pos = asarrayf( up_event_pos )
        ## The mouse must be less than 7 pixels from its down location to be considered a click.
        return mag2( down_event_pos - up_event_pos ) < 7**2
    
    def doKeyDown( self, evt ):
        ## Refresh so our status text changes.
        self.Refresh()
    def doKeyUp( self, evt ):
        ## Refresh so our status text changes.
        self.Refresh()
    
    def doDraw( self ):
        def draw_label():
            self.prepare2D()
            if self.GetMouseState().AltDown():
                draw_mode_label_text( 'Zoom out' )
            else:
                draw_mode_label_text( 'Zoom in' )
        
        def draw_rect():
            self.prepare2D()
            
            ## Draw the black selection rectangle.
            one = asarray( self.event2world( self.down_event_loc ) )
            two = asarray( self.event2world( self.GetMousePosition() ) )
            
            rect = [
                ( one[0], one[1] ),
                ( two[0], one[1] ),
                ( two[0], two[1] ),
                ( one[0], two[1] )
                ]
            
            glColor3f( 0,0,0 )
            glLineWidth( 1 )
            draw_lineloop( rect )
            
            draw_mode_label_text( 'Zoom to rect' )
        
        callback = draw_label
        ## If we're more than a click, draw the rect.
        if self.down_event_loc is not None and not self.IsClick( self.down_event_loc, self.GetMousePosition() ):
            callback = draw_rect
        
        self.DrawScene( callback = callback )

class PanMode( ModeProxy ):
    def __init__( self, canvas ):
        ModeProxy.__init__( self, canvas )
        
        self.last_event_pos = None
        
        ## TODO: Change the cursor to a hand.
        ## TODO: Preserve the selection.
    
    def doMouseDown( self, evt ):
        self.last_event_pos = tuple( evt.GetPosition() )
    
    def doMouseDrag( self, evt ):
        cur_event_pos = tuple( evt.GetPosition() )
        
        if self.last_event_pos == cur_event_pos: return
        
        ## Translate the view rect.
        self.Pan( asarray( self.event2world( self.last_event_pos ) ) - asarray( self.event2world( cur_event_pos ) ) )
        self.last_event_pos = cur_event_pos
        
        self.GetTopLevelParent().Refresh()
    
    def doDraw( self ):
        def draw_label():
            self.prepare2D()
            draw_mode_label_text( 'Drag to pan' )
        
        ## Don't draw the stroke.
        self.DrawScene( callback = draw_label )

class SameLengthAnnotation( ModeProxy ):
    def __init__( self, canvas ):
        ModeProxy.__init__( self, canvas )
    
    def doDraw( self ):
        ## Show a legend area with all annotations, plus an extra one labelled "New...".
        ## Text at bottom should say
        ## "choose / add to / remove from <color>same length group</color>",
        ## if over a spine, depending on whether
        ## no group is selected / axis is not already in the group / axis is already in the group
        ## or "start new same length group" if not over a spine.
        pass

class SweepShapeSpineDeformer( StrokeDrawer ):
    def __init__( self, canvas, shape ):
        StrokeDrawer.__init__( self, canvas, [] )
        
        assert shape is not None
        self.shape = shape
        
        self.spine_undeformed = array( self.shape.spine2d )
        self.vertex_index = None
        self.start_index = None
        self.end_index = None
    
    def doPointAdded( self ):
        ## If this is the first point, find the closest point on the spine.
        if self.vertex_index is None:
            handle = PickSweepShapeHandle(
                self.shape,
                self,
                self.stroke[0],
                parts = [SweepShapePart.Spine]
                )
            if handle is None:
                unforward_events( self )
                self.Refresh()
                return
            
            self.SaveStateForUndo()
            
            t = handle[1]
            pt = self.shape.spine2d_at_t( t )
            vertex_index = self.shape.vertex_index_at_t( t )
            self.vertex_index = vertex_index
        
        self._peel()
        #self._deform_entire_curve()
        
        self.Refresh()
    
    def _peel( self ):
        ## Peeling interface.
        ## 1 Take the difference d between the last point and the first point.
        ## 2 Find the first point along the curve in both directions farther from the
        ##   first point than d (either arc length along the curve or in world-space).
        ##   These points determine the live area of the curve.
        ## 3 If a point was found in the given direction, constrain it.
        ##   If no points were found in either direction (only one constraint at self.vertex_index),
        ##   simply translate the entire curve.
        ## 4 Deform the curve.
        
        ## 1
        d = mag( asarray( self.stroke[-1] ) - self.stroke[0] )
        
        ## 2
        start_index = self.vertex_index - 1
        while start_index >= 0:
            if line_strip_length( self.spine_undeformed[ start_index : self.vertex_index + 1 ] ) > d:
                break
            else:
                start_index -= 1
        
        end_index = self.vertex_index + 1
        while end_index < len( self.spine_undeformed ):
            if line_strip_length( self.spine_undeformed[ self.vertex_index : end_index + 1 ] ) > d:
                break
            else:
                end_index += 1
        
        #print "self.vertex_index:", self.vertex_index
        #print "start_index:", start_index
        #print "end_index:", end_index
        
        kOneEndPointsFixed = False
        kBothEndPointsFixed = True
        
        ## Make sure at least one endpoint (the one farther from self.vertex_index) stays fixed.
        if kOneEndPointsFixed:
            if len( self.spine_undeformed )-1 - self.vertex_index > self.vertex_index:
                end_index = min( len( self.spine_undeformed )-1, end_index )
            else:
                start_index = max( 0, start_index )
        
        ## Make sure both endpoints stay fixed.  (This will be overridden if self.vertex_index
        ## is an endpoint.)
        if kBothEndPointsFixed:
            end_index = min( len( self.spine_undeformed )-1, end_index )
            start_index = max( 0, start_index )
        
        ## Q: Should the live area only grow and never shrink, even if the last point moves
        ##    back closer to the first point?
        ## A: Yes, this is how ARAP and FiberMesh work.
        if self.start_index is not None: start_index = min( start_index, self.start_index )
        self.start_index = start_index
        
        if self.end_index is not None: end_index = max( end_index, self.end_index )
        self.end_index = end_index
        
        #print "self.start_index:", self.start_index
        #print "self.end_index:", self.end_index
        
        ## 3
        ## If start/end_index fall off the end, don't constrain them.
        constrain_start_index = start_index >= 0
        tangent_constrain_start_index = start_index - 1 >= 0
        start_index = max( start_index - 1, 0 )
        
        constrain_end_index = end_index < len( self.spine_undeformed )
        tangent_constrain_end_index = end_index + 1 < len( self.spine_undeformed )
        end_index = min( end_index + 1, len( self.spine_undeformed ) - 1 )
        
        constraints = {}
        if constrain_start_index:
            constraints[ 0 ] = self.spine_undeformed[ start_index ]
        if tangent_constrain_start_index:
            constraints[ 1 ] = self.spine_undeformed[ start_index + 1 ]
        if constrain_end_index:
            constraints[ end_index - start_index ] = self.spine_undeformed[ end_index ]
        if tangent_constrain_end_index:
            constraints[ end_index - 1 - start_index ] = self.spine_undeformed[ end_index - 1 ]
        
        constraints[ self.vertex_index - start_index ] = self.stroke[ -1 ]
        
        if len( constraints ) == 1:
            delta = ( asarray( self.stroke[ -1 ] ) - self.stroke[ 0 ] )
            self.shape.set_spine2d( self.spine_undeformed + delta )
            return
        
        ## 4
        from laplacian_curve_editing import new_curve_from_curve_and_constraints
        new_curve = array( self.spine_undeformed )
        new_curve[ start_index : end_index + 1 ] = new_curve_from_curve_and_constraints(
            self.spine_undeformed[ start_index : end_index + 1 ],
            curve_is_closed = False,
            constraints = list(constraints.items())
            )
        self.shape.set_spine2d( new_curve )
    
    def _deform_entire_curve( self ):
        constraints = {}
        constraints[ 0 ] = self.spine_undeformed[ 0 ]
        constraints[ len( self.spine_undeformed ) - 1 ] = self.spine_undeformed[ -1 ]
        
        constraints[ self.vertex_index ] = self.stroke[ -1 ]
        
        from laplacian_curve_editing import new_curve_from_curve_and_constraints
        self.shape.set_spine2d( new_curve_from_curve_and_constraints(
            self.spine_undeformed,
            curve_is_closed = False,
            constraints = list(constraints.items())
            ) )
    
    def doStrokeFinished( self ):
        self.shape.smooth_spine2d( 0 )
        self.annotations.Beautify()
        self.GetTopLevelParent().Refresh()
        
        self.shape = None
    
    def doMouseMove( self, evt ):
        StrokeDrawer.doMouseMove( self, evt )
        
        ## Refresh so our Smart Guide changes.
        self.Refresh()
    
    def doDraw( self ):
        # print( "SweepShapeSpineDeformer.doDraw()" )
        def draw_callback():
            # print( "SweepShapeSpineDeformer.doDraw().draw_callback()" )
            self.prepare2D()
            
            ## Draw the live area of spine in red.
            if self.start_index is not None and self.end_index is not None:
                self.glLineWidthWorld( WorldDefaultLineWidth( self ) )
                glColor3f( 1, 0, 0 )
                draw_linestrip( self.shape.spine2d[ max( self.start_index, 0 ) : min( self.end_index + 1, len( self.shape.spine2d ) ) ] )
            
            else:
                ## Draw a Smart Guide
                assert len( self.stroke ) == 0
                DrawPickedSweepShapeHandle(
                    self.shape,
                    PickSweepShapeHandle(
                        self.shape,
                        self,
                        self.event2world( self.GetMousePosition() ),
                        parts = [SweepShapePart.Spine]
                        )
                    )
            
            draw_mode_label_text( 'Drag to deform' )
        
        self.DrawScene( callback = draw_callback )
