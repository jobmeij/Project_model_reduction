function T = CalculateOutput(model,phi,sys)
T = phi.x*sys.a;
T = T + model.Tamb;
end