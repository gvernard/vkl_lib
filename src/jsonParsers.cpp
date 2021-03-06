#include "jsonParsers.hpp"
#include "lightProfile.hpp"
#include "massModels.hpp"

#include <map>
#include <iostream>
CollectionMassModels JsonParsers::parse_mass_model(const Json::Value json,std::string prefix){
  CollectionMassModels mycollection;
  mycollection.models.resize(json.size());
  for(int k=0;k<json.size();k++){
    std::string type = json[k]["type"].asString();
    std::map<std::string,std::string> pars;
    for(Json::Value::const_iterator it=json[k]["pars"].begin();it!=json[k]["pars"].end();it++){
      pars.insert( std::pair<std::string,std::string> (it.key().asString(),it->asString()) );      
    }
    if( type == "pert" ){
      pars["filepath"] = prefix + pars["filepath"];
    }
    mycollection.models[k] = FactoryMassModel::getInstance()->createMassModel(type,pars);
  }
  return mycollection;
}
  
CollectionProfiles JsonParsers::parse_profile(const Json::Value json,std::string prefix){
  CollectionProfiles mycollection;
  mycollection.profiles.resize(json.size());
  for(int k=0;k<json.size();k++){
    std::string type = json[k]["type"].asString();
    std::map<std::string,std::string> pars;
    for(Json::Value::const_iterator it=json[k]["pars"].begin();it!=json[k]["pars"].end();it++){
      pars.insert( std::pair<std::string,std::string> (it.key().asString(),it->asString()) );      
    }
    if( type == "custom" ){
      pars["filepath"] = prefix + pars["filepath"];
    }
    mycollection.profiles[k] = FactoryProfile::getInstance()->createProfile(type,pars);
  }
  return mycollection;
}
