void DrawMicromegas_tgeom()
{
    
    const TString inputFile = "sPHENIXGeom.root";
    
    // open TFile
    auto f = TFile::Open( inputFile );
    
    // get geomanager
    auto geomanager = static_cast<TGeoManager*>( f->Get("Default" ) );

    auto volumes = geomanager->GetListOfVolumes();
    for( const auto&& object:*volumes )
    {
        // cast to volume
        const auto volume = static_cast<TGeoVolume*>( object );
        
        // get name
        const TString name( volume->GetName() );
        // const bool visible = name.BeginsWith( "MICROMEGAS" );
        const bool visible = name.BeginsWith( "MICROMEGAS" ) || name.BeginsWith( "tpc" );
        volume->SetVisibility( visible );
    }
    
    // draw master volume
    geomanager->GetTopVolume()->Draw( "ogl" );
}
