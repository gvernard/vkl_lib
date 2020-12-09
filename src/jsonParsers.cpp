#include "jsonParsers.hpp"
#include "sourceProfile.hpp"
#include "massModels.hpp"

#include <map>

CollectionMassModels JsonParsers::parse_mass_model(const Json::Value json){
  CollectionMassModels mycollection;
  mycollection.models.resize(json.size());
  for(int k=0;k<json.size();k++){
    std::string type = json[k]["type"].asString();
    std::map<std::string,std::string> pars;
    for(Json::Value::const_iterator it=json[k]["pars"].begin();it!=json[k]["pars"].end();it++){
      pars.insert( std::pair<std::string,std::string> (it.key().asString(),it->asString()) );      
    }
    mycollection.models[k] = FactoryMassModel::getInstance()->createMassModel(type,pars);
  }
  return mycollection;
}
  
CollectionProfiles JsonParsers::parse_profile(const Json::Value json){
  CollectionProfiles mycollection;
  mycollection.profiles.resize(json.size());
  for(int k=0;k<json.size();k++){
    std::string type = json[k]["type"].asString();
    std::map<std::string,std::string> pars;
    for(Json::Value::const_iterator it=json[k]["pars"].begin();it!=json[k]["pars"].end();it++){
      pars.insert( std::pair<std::string,std::string> (it.key().asString(),it->asString()) );      
    }
    mycollection.profiles[k] = FactoryProfile::getInstance()->createProfile(type,pars);
  }
  return mycollection;
}
