// Version menu for the published documentation. The site has one directory
// per version (latest, v0.0.2, ...) and a versions.json listing them; a
// local build without that file shows no menu.
(function () {
  "use strict";
  function current() {
    var parts = window.location.pathname.split("/");
    return parts.length > 2 ? parts[parts.length - 2] : "";
  }
  function build(versions) {
    var here = current();
    if (versions.indexOf(here) < 0) return;
    var select = document.createElement("select");
    select.setAttribute("aria-label", "Documentation version");
    versions.forEach(function (name) {
      var option = document.createElement("option");
      option.value = name;
      option.textContent = name;
      option.selected = name === here;
      select.appendChild(option);
    });
    select.addEventListener("change", function () {
      window.location.href = "../" + select.value + "/index.html";
    });
    var box = document.createElement("div");
    box.style.margin = "0.5em 0";
    box.appendChild(document.createTextNode("Version: "));
    box.appendChild(select);
    var anchor = document.querySelector(".sphinxsidebarwrapper") || document.body;
    anchor.insertBefore(box, anchor.firstChild);
  }
  document.addEventListener("DOMContentLoaded", function () {
    fetch("../versions.json")
      .then(function (response) { return response.ok ? response.json() : []; })
      .then(function (versions) { if (Array.isArray(versions)) build(versions); })
      .catch(function () {});
  });
})();
