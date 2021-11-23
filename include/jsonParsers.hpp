#ifndef JSON_PARSERS_HPP
#define JSON_PARSERS_HPP

#include "json/json.h"

class CollectionMassModels;
class CollectionProfiles;

class JsonParsers {
public:
  static CollectionMassModels parse_mass_model(const Json::Value mass_model,std::string prefix="");
  static CollectionProfiles parse_profile(const Json::Value profile,std::string prefix="");
  static CollectionProfiles parse_profile(const Json::Value profile,double ZP,std::string prefix="");
};

#endif /* JSON_PARSERS_HPP */
