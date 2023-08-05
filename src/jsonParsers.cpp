#include "jsonParsers.hpp"
#include "lightProfile.hpp"
#include "massModels.hpp"

#include <map>
#include <iostream>

using namespace vkl;

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
    if( json[k].isMember("mass-to-light") ){
      for(Json::Value::const_iterator it=json[k]["mass-to-light"].begin();it!=json[k]["mass-to-light"].end();it++){
	pars.insert( std::pair<std::string,std::string> (it.key().asString(),it->asString()) );
      }
    }
    mycollection.profiles[k] = FactoryProfile::getInstance()->createProfile(type,pars);
  }
  return mycollection;
}

CollectionProfiles JsonParsers::parse_profile(const Json::Value json,double ZP,std::string prefix){
  // Same as without the ZP. Here I just add the ZP to the profile pars.
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
    if( json[k].isMember("mass-to-light") ){
      for(Json::Value::const_iterator it=json[k]["mass-to-light"].begin();it!=json[k]["mass-to-light"].end();it++){
	pars.insert( std::pair<std::string,std::string> (it.key().asString(),it->asString()) );
      }
    }
    pars.insert( std::pair<std::string,std::string> ("ZP",std::to_string(ZP)) ); // Adding the ZP
    mycollection.profiles[k] = FactoryProfile::getInstance()->createProfile(type,pars);
  }
  return mycollection;
}
