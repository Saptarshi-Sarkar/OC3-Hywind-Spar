function Platform = CreatePlatform()
% Initiate Hydrodynamic Load Parameters
TotDraft             = 120;
TaprTop              = 4;
TaprBot              = 12;
TaprTopDia           = 6.5;
TaprBotDia           = 9.4;
Platform.nStrips     = 240;
Platform.StripWidth  = TotDraft/Platform.nStrips;
Platform.StripDepths = linspace(0,TotDraft,Platform.nStrips);
Platform.StripDia    = interp1([0 TaprTop TaprBot TotDraft], [TaprTopDia TaprTopDia TaprBotDia TaprBotDia], Platform.StripDepths);
Platform.AM          = AddedMass( Platform );
